"""

BurriTE: A wrapper pipeline to annotate and analyze TEs in pangenomic context V.0.4.0
Novelties (changes since v.0.0.0):

  === GraffiTE container/runtime handling (v.0.4.0) ===
  * --graffite_image optionally supplies a previously downloaded local SIF/IMG;
    omitting it preserves GraffiTE's configured default container behaviour.
  * --graffite_tmpdir optionally supplies a shared writable temporary directory for
    GraffiTE tasks. Its default is <output>/01_GraffiTE/tmp, avoiding ephemeral SLURM
    TMPDIR paths that may not exist inside the container or on other compute nodes.
  * TMP/TMPDIR are injected into the GraffiTE container through APPTAINERENV_ variables
    in an isolated subprocess environment; stale SINGULARITYENV_TMP* values are removed,
    and the shared temporary directory is explicitly bind-mounted into the container.
  * GraffiTE is launched with -resume, and stdout/stderr are always written to
    01_GraffiTE/GraffiTE_log.txt. Non-zero Nextflow exit codes now stop BurriTE early.

  === Output directory layout (numbered by step order) ===
  * Intermediate/result folders are numbered in the order they are produced:
      01_GraffiTE  (SV discovery)         02_RepeatMasker (RM + OneCode per genome)
      03_annotation_gff (RM/GraffiTE gff) 04_graffite_liftover (GraffiTE ref->ALT)
      05_joined_annotation (per-genome RM+GraffiTE = <sp>_final_annot.gff)
      06_lift_to_ref (filter/flanks/minimap + the FIRST lift, ALT anchors -> reference)
      07_merge (step 11 consensus collapse + step 12/13 flank liftback intermediates)
      08_liftback (step 14 <g>_liffout_polished.gff + <g>_Transfers_TE.txt genotype tables)
      09_final (deliverables: Presence_Absence_matrix.txt + per_genome_gff/)

  === Merge / collapse (step 11) ===
  * The per-genome-to-reference merge that failed in the original at the bedtools
    groupby step now works. It uses the Yoss's run_pipeline.sh logic and requires
    bedtools >= 2.31 (2.26 silently collapses every row of a group onto one line;
    the original was written against a broken groupby). Run inside an env where
    `bedtools --version` reports 2.31.x.
  * --merge_method {groupby, incremental} chooses how per-genome TE calls are merged
    into a non-redundant consensus catalogue; incremental is now the GLOBAL DEFAULT:
      - incremental (default): Python port of yoss_scripts/incremental_collapse/
        build_incremental_collapse.sh. Adds genomes one at a time and merges an
        interval into an existing cluster only if same chromosome + same key (see
        --merge_by) + coordinate overlap + cluster internal length spread stays within
        --lenthr. Keeps length-variant and private loci that the groupby merges away,
        labels each consensus locus shared (>=2 genomes) or singleton (1), gives loci
        opaque ids LOC<n>, and records a Nested flag (TRUE if a locus is fully contained
        in another consensus locus of a different key = true TE-in-TE nesting).
      - groupby (opt-in legacy): superfamily-level collapse (original behaviour, repaired).
  * --merge_by {family, superfamily} (default superfamily; incremental only) sets the
    clustering granularity: superfamily = Class-level (e.g. DNA_RC); family = Name-level
    (e.g. ATREP13, stricter). family and superfamily are written as separate matrix columns.
  * CLI knobs for the incremental method: --lenthr (100), --minlen (0), --bedmin (100).
    minlen drops absence-derived degenerate input intervals; bedmin drops sub-100 bp
    consensus loci from the second-mapping BED so phase 4 does not genotype 1 bp "TEs".

  === TE-copy length filtering (--min_length / --max_length) ===
  * The step-7 annotation filter (filter_te_table, plus build_reference_bed and the
    --final_annotation all overlay) keeps a TE copy only if min_length <= len <= max_length.
    Both are now UNRESTRICTED by default: --min_length default 0, --max_length default unset
    (no upper cap, handled internally as "no cap"). Previously the hard 100/15000 bp window
    silently discarded very short and very long copies (long LTR elements were clipped at
    15 kb); with no cap those copies flow through discovery, collapse and genotyping.

  === Real TE length carried through (fix to BOTH methods) ===
  * The per-species BED now carries LengthT (real TE size, OneCode col 7) instead
    of LengthR (the mapped flank gap, col 8): the awk emits `$6"|"$7`, not `...|$8`.
    Previously the consensus label length was the distance between the translated
    flanks, which for absence-derived / rearranged loci degenerates to a few bp and
    produced spurious 7 bp "TEs".
  * The collapse/cluster labelling keys on LengthT. Consensus loci are labelled
    `FM<idx>_<superfamily>|<median LengthT>`; the old coordinate-span drop filter was
    removed so length-consistent groups are kept.

  === Phase 4 transfer-back genotyping (steps 12-14) repaired ===
    (these were pre-existing-broken in BurriTE_original and now run end to end)
  * Step 12 generate_flanking_gffs now reads the 4-column v2 consensus BED
    (chr,start,end,name) and derives the label from name.split("|")[0], instead of
    assuming a 6-column input (which raised a float64 TypeError on `label + "_F1"`).
  * Step 13 lifts the per-chromosome consensus flanks (07_merge/FM_{chrm}.gff
    with {chrm}_gff_features.txt) back onto every assembly; output prefix
    07_merge/{assembly}_{chrm}. PATH now includes the burrite conda bin so
    `liftoff` is found.
  * Step 14 concatenates the per-chromosome polished GFFs per assembly and guards
    against missing per-chromosome files (the original crashed on `(f1 or f2)[0]`
    when both flanks were absent; now `(f1 or f2 or ("NA",))[0]`).
  * evaluate_te_transfer implements the validated Present/Absent definition:
      gap = F2_start - F1_end
      Absent    if 0 <= gap <= 10
      Present   if |LengthT - gap| <= 0.20*LengthT  (tolerance scales with TE length; LengthT read from the consensus BED label)
      Ambiguous otherwise
      NA        if a flank did not lift
    and emits an Unmappable_flanks indicator column. The tolerance is now RELATIVE
    (0.20*LengthT) instead of an absolute +/-50 bp; the same relative rule gates the
    step-10 first-transfer length filter and the reference self-genotype, so long TEs
    are no longer forced into Ambiguous just because a fixed 50 bp window is a tiny
    fraction of their length.

  === Reference TE annotation integrated (opt-in) ===
  * New --include_reference {Y,N} (no default). The reference's own joined
    RepeatMasker+GraffiTE annotation was previously computed but never consumed, so
    reference(-private) TEs never entered the catalogue and were never genotyped. When
    Y, build_reference_bed reformats that annotation (slot 7) into a per-genome BED in
    bed_files/ (reference TEs are already in reference coordinates, so no ALT->reference
    liftoff is needed) and it is merged/collapsed at step 11 together with the ALT-lifted
    loci; overlapping reference copies collapse against the ALT loci, and genuinely
    reference-extra loci become new consensus rows that steps 12-15 then genotype in
    every query. --merge_method defaults to incremental in this mode (precise
    family+overlap+length clustering of the many native reference copies).
  * evaluate_reference_transfer genotypes the reference itself at each consensus locus
    (gap = FM end - start, the identity liftback) and step 15 adds it as the first
    matrix column, so the matrix is no longer reference-blind.

  === Annotation-Present override (always on) ===
  * The step-11 collapse already records which genomes' OWN annotation contributed each
    consensus locus (incremental: annots_incremental_collapsed.tsv samples column; groupby:
    annots_groupby.bed col 6). build_contributor_map exposes this for BOTH merge methods,
    and step 15 marks a contributing genome `Present:annot` at that locus whenever its
    flank-liftback genotype was NOT already Present (Absent/Ambiguous/NA). Rationale: direct
    annotation is more reliable than the query->ref->query flank round-trip, which
    mis-genotypes short/repetitive TEs; this rescues those and removes all "present in 0
    genomes" loci (every consensus locus has >=1 contributing genome). Present cells are
    untouched and the override is labelled (auditable, like NA:<reason>); anything
    startswith "Present" is present.

  === New outputs ===
  * Step 15 presence/absence matrix (09_final/Presence_Absence_matrix.txt): one row per
    consensus locus, columns `Chrm start end LOC_ID family superfamily TE_length <genomes>`;
    cells Present/Absent/Ambiguous/NA (NA annotated `NA:<reason>`), plus Present:annot from
    the annotation-Present override.
  * Step 16 per-genome GFFs (09_final/per_genome_gff/), controlled by --final_annotation:
      - burrite (default): <g>_TEs.gff = the consensus TEs Present or Absent in genome g
        (Ambiguous and Not_mapped are excluded). Feature type encodes the focal state
        (transposable_element = present copy, insertion_site = absent site); attributes carry
        ID=LOC_ID, family, superfamily, coord_source, and four tags (Present/Absent/Ambiguous/
        Not_mapped) listing the OTHER genomes in each state ("none" if empty).
      - all: <g>_all_TEs.gff additionally folds in genome g's own RepeatMasker + GraffiTE
        annotations (length-filtered like step 7, remapped to reference-chromosome coords),
        de-duplicated against the BurriTE loci and each other by reciprocal overlap (>= the
        --final_annotation_dedup_overlap fraction, default 0.5) with
        precedence BurriTE > GraffiTE > RepeatMasker; the surviving feature records
        merged_sources. RM/GraffiTE features get the four state tags = not_applicable.
    Coordinates (coord_source attribute): confident genotypes use the real inter-flank span
    ("flank_gap"); Present:annot overrides take the EXACT own-genome annotation coordinates of
    the copy that genome contributed, recovered by joining the step-11 collapse dictionary to
    that genome's own-coordinate annotation bed ("annotation"), which resolves the paralog
    ambiguity and avoids the unreliable flank-liftback gap (fallback "annot_anchor" = a lifted
    flank + the known TE length, used only when the dictionary lookup is unavailable).

  === Robustness fixes ===
  * te_transfer_per_species uses os.path.abspath() on the input GFF and output dir,
    so the internal `cd` no longer breaks relative gff paths (original failed with
    `grep: ..._liffout.gff_polished: No such file`).
  * Species tag is now assembly.replace("_","") before it is embedded in TE IDs, so
    genomes like C_arabica no longer split the superfamily/species fields during the
    merge (which previously stopped cross-genome loci from collapsing).
  * parse_liftoff_polished prefers copy_num_ID= then falls back to ID= when reading
    lifted feature identifiers.

  === Parallel liftoff (steps 10 and 13) ===
  * The per-(assembly x chromosome) liftoff calls in transfer_anchors (step 10) and the
    step-13 liftback previously ran one after another. They are now dispatched to a thread
    pool of LIFTOFF_MAX_PARALLEL (=10) concurrent jobs; each liftoff still gets its own
    internal -p threads (cores // LIFTOFF_MAX_PARALLEL).
  * Concurrency-safety: liftoff/minimap2/pysam write sidecar index files (.fai, .mmi) next
    to the input fasta paths, so jobs that share a reference fasta would race on them --
    corrupting one job's index (truncated .fai = crash; clobbered .mmi = a silently EMPTY
    lift). Two guards prevent this: (a) ensure_faidx pre-builds the .fai serially before the
    pool, and (b) run_liftoff points each job at PRIVATE per-job symlinks (target.fa/
    reference.fa) inside its own unique -dir, so every sidecar index is unique to the job.

"""

import sys
import os
import io
import re
import argparse
import pandas as pd
import multiprocessing
from Bio import SeqIO
import subprocess
try:
    import psutil
except ImportError:  # Stage helpers do not require it; the full Conda env includes it.
    psutil = None
import shutil
from pathlib import Path
from collections import defaultdict
import concurrent.futures

# Step 10 (transfer_anchors) and step 13 run one liftoff per (assembly x chromosome). These are
# independent, so they are dispatched to a thread pool of this many concurrent liftoff jobs instead
# of running sequentially. Each liftoff still uses its own -p (internal) parallelism; see run_liftoff.
LIFTOFF_MAX_PARALLEL = 10

def ensure_faidx(fasta):
    """Build the samtools .fai index for a fasta SERIALLY (up front) if missing. Concurrent liftoff
    jobs that share a reference fasta would otherwise each try to build the .fai at the same time and
    race, corrupting it ('OSError: truncated file'). Pre-building makes the concurrent jobs read-only."""
    if fasta and os.path.exists(fasta) and not os.path.exists(fasta + ".fai"):
        subprocess.run(["samtools", "faidx", fasta],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def write_sequences_file(sequences, filename):
    try:
        SeqIO.write(sequences, filename, "fasta")
    except FileNotFoundError:
        print("FATAL ERROR: I couldn't find the file, please check: '" + filename + "'. Path not found")
        sys.exit(0)
    except PermissionError:
        print("FATAL ERROR: I couldn't access the files, please check: '" + filename + "'. I don't have permissions.")
        sys.exit(0)
    except Exception as exp:
        print("FATAL ERROR: There is a unknown problem writing sequences in : '" + filename + "'.")
        print(exp)
        sys.exit(0)


def create_output_folders(folder_path):
    try:
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
    except FileNotFoundError:
        print("FATAL ERROR: I couldn't create the folder " + folder_path + ". Path not found")
        sys.exit(0)
    except PermissionError:
        print("FATAL ERROR: I couldn't create the folder " + folder_path + ". I don't have permissions.")
        sys.exit(0)


def delete_files(file_path):
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except PermissionError:
            print("WARNING: The file " + file_path + " couldn't be removed because I don't have permissions.")


def run_GraffiTE(assemblies, TE_lib, reference_assembly, cores, outputdir, verbose,
                 graffite_image=None, graffite_tmpdir=None):
    create_output_folders(outputdir)

    if not os.path.exists(outputdir + "/3_TSD_search/pangenome.vcf"):
        if graffite_tmpdir is None:
            graffite_tmpdir = os.path.join(outputdir, "tmp")
        graffite_tmpdir = os.path.abspath(os.path.expanduser(graffite_tmpdir))

        try:
            os.makedirs(graffite_tmpdir, mode=0o700, exist_ok=True)
        except OSError as exp:
            raise RuntimeError(
                "Could not create the GraffiTE temporary directory: "
                + graffite_tmpdir
            ) from exp

        if not os.path.isdir(graffite_tmpdir):
            raise RuntimeError(
                "The GraffiTE temporary path is not a directory: "
                + graffite_tmpdir
            )
        if not os.access(graffite_tmpdir, os.W_OK | os.X_OK):
            raise PermissionError(
                "The GraffiTE temporary directory is not writable: "
                + graffite_tmpdir
            )

        # Use a private environment for GraffiTE so that the remaining BurriTE
        # steps retain the user's original TMP/TMPDIR configuration.
        graffite_env = os.environ.copy()
        graffite_env.pop("SINGULARITYENV_TMP", None)
        graffite_env.pop("SINGULARITYENV_TMPDIR", None)
        graffite_env.update({
            "TMP": graffite_tmpdir,
            "TMPDIR": graffite_tmpdir,
            "APPTAINERENV_TMP": graffite_tmpdir,
            "APPTAINERENV_TMPDIR": graffite_tmpdir,
        })

        # The temporary path may live outside Nextflow's task work directory.
        # Bind it explicitly so it exists at the same absolute path in the
        # Apptainer container, while preserving any user-defined bind paths.
        existing_bindpath = graffite_env.get("APPTAINER_BINDPATH", "")
        graffite_env["APPTAINER_BINDPATH"] = (
            existing_bindpath + "," + graffite_tmpdir
            if existing_bindpath else graffite_tmpdir
        )

        print("MESSAGE: GraffiTE temporary directory: " + graffite_tmpdir)

        command = [
            'nextflow', 'run', 'cgroza/GraffiTE',
            '-profile', 'cluster', '-resume'
        ]
        if graffite_image is not None:
            # An absolute path is required so the image remains resolvable from
            # Nextflow work directories and from the cluster compute nodes.
            command.extend(['-with-singularity', graffite_image])
        command.extend([
            '--assemblies', assemblies, '--TE_library', TE_lib,
            '--reference', reference_assembly, '--graph_method', 'graphaligner',
            '--genotype', 'false', '--cores', str(cores), '--out', outputdir
        ])

        output = subprocess.run(
            command,
            env=graffite_env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if verbose:
            log_path = outputdir + "/GraffiTE_log.txt"
            with open(log_path, "w") as logfile:
                logfile.write(output.stdout)
                logfile.write("\n\nSTDERR:\n")
                logfile.write(output.stderr)

            if output.returncode != 0:
                raise RuntimeError(
                    "GraffiTE failed with exit status " + str(output.returncode)
                    + ". See: " + log_path
                )

    else:
        print("WARNING: GraffiTE was already run. Skipping....")

    return outputdir + "/3_TSD_search/pangenome.vcf"


def run_RepeatMasker(file_dict, TE_lib, outputdir, cores, verbose):

    for assembly in file_dict.keys():
        create_output_folders(outputdir + '/RM_'+assembly)
        nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][0]))
        if not os.path.exists(outputdir + '/RM_'+assembly+'/'+''.join(nombre_archivo)+'.out'):
            output = subprocess.run(
                ['RepeatMasker', '-famdb_dir', '', '-pa', str(cores), '-dir', outputdir + '/RM_'+assembly, '-lib', TE_lib, file_dict[assembly][0]],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            file_dict[assembly].append(outputdir + '/RM_'+assembly+'/'+''.join(nombre_archivo)+'.out')

            if verbose:
                with open(outputdir + "/RepeatMasker_log.txt", "w") as logfile:
                    logfile.write(output.stdout)
                    logfile.write("\n\nSTDERR:\n")
                    logfile.write(output.stderr)
        else:
            print("WARNING: RepeatMasker already found in the path: " + outputdir + '/RM_'+assembly+'/'+''.join(nombre_archivo)+'.out. Skipping....')
            file_dict[assembly].append(outputdir + '/RM_' + assembly + '/' + ''.join(nombre_archivo) + '.out')
    return file_dict


def run_OneCode(file_dict, outputdir, tools_path, verbose):
    for assembly in file_dict.keys():
        RM_output = file_dict[assembly][1]

        if not os.path.exists(RM_output+'.elem_sorted.csv.out'):
            output1 = subprocess.run(
                [tools_path+'/OneCodeToFindThemAll/build_dictionary.pl', '--rm',
                 RM_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            with open(outputdir + '/RM_'+assembly+'/dico.txt', "w") as dico:
                dico.write(output1.stdout)
                logfile.write("\n\nSTDERR:\n")
                logfile.write(output1.stderr)

            output2 = subprocess.run(
                [tools_path+'/OneCodeToFindThemAll/one_code_to_find_them_all_but_sanely.pl', '--rm',
                RM_output, '--ltr', outputdir + '/RM_' + assembly + '/dico.txt', '--unknown',
                 '--strict'], stdout = subprocess.PIPE, stderr=subprocess.PIPE, text = True)

            if verbose:
                with open(outputdir + '/RM_'+assembly+"/OneCode_log.txt", "w") as logfile:
                    logfile.write("####### Command 2 #######\n\n" + output2.stdout)
                    logfile.write("\n\nSTDERR:\n")
                    logfile.write(output2.stderr)

            annotations = []
            with open(RM_output+".elem_sorted.csv", "r") as oneCode:
                lines = oneCode.readlines()
                for line in lines:
                    if line.startswith("###"):
                        annotations.append(line.replace("###", ""))

            with open(RM_output + ".elem_sorted.csv.out", "w") as oneCode:
                oneCode.write("".join(annotations))

            file_dict[assembly].append(RM_output+'.elem_sorted.csv.out')
        else:
            print("WARNING: OneCode already found in the path: " + outputdir + '/RM_'+assembly + ". Skipping ...")
            file_dict[assembly].append(RM_output + '.elem_sorted.csv.out')
    return file_dict


def parse_graffite_to_repeatmasker(file_dict, reference, input_file, outputdir):
    species = {}
    with open(input_file, "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue  # Omitir lÃ­neas incompletas

                species_name = fields[2].split(".")
                SV_type = species_name[-2]

                if SV_type == "DEL":
                    # a deletion means a TE in Ref not present in the ALT assemblies
                    species_name = reference
                elif SV_type == "INS":
                    # an insertion means a TE not present in Ref but present in the ALT assemblies
                    species_name = species_name[0]  # + "." + species_name[1]

                contig = fields[0]
                start = int(fields[1])
                info = fields[7]
                genotype = fields[9] if len(fields) > 9 else "."

                # Extraer informaciÃ³n relevante del campo INFO
                end_match = re.search(r"SVLEN=(-*\d+)", info)
                repeat_match = re.search(r"repeat_ids=([\w\-_]+)", info)
                class_match = re.search(r"matching_classes=([\w/\-]+)", info)
                strand_match = re.search(r"STRANDS=([+-])", info)
                end = int(end_match.group(1)) if int(end_match.group(1)) > 0 else int(end_match.group(1))*-1
                end = start + end
                repeat_name = repeat_match.group(1) if repeat_match else "Unknown"
                repeat_class = class_match.group(1) if class_match else "Unknown"
                strand = strand_match.group(1) if strand_match else "+"

                # Formato RepeatMasker (.out)
                sw_score = SV_type  # Sin informaciÃ³n, se puede modificar
                divergence, deletions, insertions = 0, 0, 0  # Valores por defecto
                output_line = (
                    f"{sw_score}\t{divergence}\t{deletions}\t{insertions}\t{contig}\t{start}\t{end}\t{end - start}\t{strand}\t{repeat_name}\t{repeat_class}\t"
                    f"{0}\t{0}\t{0}\t{0}\t{0}\t{'No_ref_available'}\n")

                if species_name in species.keys():
                    species[species_name].append(output_line)
                else:
                    species[species_name] = [output_line]

    for species_name in species:
        with open(outputdir + "/" + species_name + "_GraffiTE.out", "w") as outfile:
            outfile.write("\n".join(species[species_name]))
        file_dict[species_name].append(outputdir + "/" + species_name + "_GraffiTE.out")

    return file_dict


def RM_out_2_gff(file_dict, outputdir):
    for assembly in file_dict.keys():
        for i in range(2, 4, 1):  # to do both: RepeatMasker (i=2) and GraffiTE.out (i=3)
            input_file = file_dict[assembly][i]
            nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][i]))[0]

            # Cargar el archivo tabulado (.out) generado por RepeatMasker (sin cabeceras de texto)
            df = pd.read_csv(input_file, sep="\t", header=None, comment="#")

            # Asumimos que las columnas son estas (puedes ajustar si el orden cambia)
            df.columns = [
                "SW_score", "perc_div", "perc_del", "perc_ins",
                "query_seq", "query_start", "query_end", "length",
                "sense", "element", "family", "Pos_Repeat_Beg", "Pos_Repeat_End",
                "Pos_Repeat_Left", "ID", "Num_Assembled", "%_of_Ref"
            ]

            # Asegurar que las posiciones estÃ¡n en el orden correcto
            df["start"] = df[["query_start", "query_end"]].min(axis=1)
            df["end"] = df[["query_start", "query_end"]].max(axis=1)

            # Normalizar strand
            df["sense"] = df["sense"].replace({"C": "-"})  # "C" en RepeatMasker significa complemento
            df["sense"] = df["sense"].fillna("+")          # Suponemos que si estÃ¡ vacÃ­o es "+"

            # Generar la columna de atributos para GFF3
            df["attributes"] = (
                "ID=TE" + df.index.astype(str) + ";" +
                "Name=" + df["element"] + ";" +
                "Class=" + df["family"]
            )

            # Crear el DataFrame GFF3
            gff3_df = pd.DataFrame({
                "seqid": df["query_seq"],
                "source": "RepeatMasker",
                "type": "repeat_region",
                "start": df["start"],
                "end": df["end"],
                "score": ".",
                "strand": df["sense"],
                "phase": ".",
                "attributes": df["attributes"]
            })

            # Guardar en GFF3
            with open(outputdir+"/"+nombre_archivo+".gff", "w") as f:
                gff3_df.to_csv(f, sep="\t", header=False, index=False)

            file_dict[assembly].append(outputdir+"/"+nombre_archivo+".gff")

    return file_dict


def transfer_graffite_to_assemblies(file_dict, reference, outputdir, verbose):
    create_output_folders(outputdir)
    reference_assembly = file_dict[reference][0]

    # To create the features file needed by LiftOff
    with open(outputdir+"/features_file.txt", "w") as fea_file:
        fea_file.write("repeat_region")

    for assembly in file_dict.keys():
        GraffiTE_RM_output = file_dict[assembly][5]
        assembly_fasta = file_dict[assembly][0]
        nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][5]))[0]
        if assembly != reference:
            if not os.path.exists(outputdir+"/"+nombre_archivo+"_lifted.gff"):
                output = subprocess.run(
                    ["liftoff", "-g", GraffiTE_RM_output, "-o", outputdir+"/"+nombre_archivo+"_lifted.gff",
                          "-f", outputdir+"/features_file.txt", assembly_fasta, reference_assembly], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

                file_dict[assembly].append(outputdir+"/"+nombre_archivo+"_lifted.gff")

                if verbose:
                    with open(outputdir + "/"+assembly+"_liftoff_log.txt", "w") as logfile:
                        logfile.write(output.stdout)
                        logfile.write("\n\nSTDERR:\n")
                        logfile.write(output.stderr)

            else:
                print("WARNING: Lifted file already found in the path: " + outputdir+"/"+nombre_archivo+"_lifted.gff. Skipping ...")
                file_dict[assembly].append(outputdir + "/" + nombre_archivo + "_lifted.gff")
        else:
            # Copying the GraffiTE GFF file for the reference genome
            shutil.copy(GraffiTE_RM_output, outputdir+"/"+nombre_archivo+"_lifted.gff")
            file_dict[assembly].append(outputdir + "/" + nombre_archivo + "_lifted.gff")

    return file_dict


def join_TE_annotations(file_dict, outputdir):
    create_output_folders(outputdir)
    for species in file_dict.keys():

        if not os.path.exists(outputdir+"/"+species+"_final_annot.gff"):
            RepeatMasker = pd.read_csv(file_dict[species][4], sep='\t', header=None)
            GraffiTE = pd.read_csv(file_dict[species][6], sep='\t', comment="#", header=None)
            joined_annot = pd.concat([RepeatMasker, GraffiTE])
            joined_annot.columns = ["seqid", "source", "type", "start",
                "end", "score", "strand", "phase", "attributes"]
            joined_annot = joined_annot.sort_values(by=['seqid', 'start', 'score'], ascending=[True, True, True])

            # Lista para almacenar anotaciones no solapadas
            merged_annotations = []
            prev_chrom, prev_start, prev_end, prev_row = None, None, None, None

            for _, row in joined_annot.iterrows():
                chrom, start, end = row['seqid'], row['start'], row['end']

                # Si no hay solapamiento, agregar la anotaciÃ³n
                if prev_chrom != chrom or start > prev_end:
                    merged_annotations.append(row)
                    prev_chrom, prev_start, prev_end, prev_row = chrom, start, end, row
                else:
                    # Si hay solapamiento, conservar la mejor anotaciÃ³n (ya estÃ¡ ordenado)
                    continue

            joined_annot.to_csv(outputdir+"/"+species+"_final_annot.gff", header=False, sep='\t', index=False)
            file_dict[species].append(outputdir+"/"+species+"_final_annot.gff")
        else:
            print("WARNING: Joined annotations already found in the path: " + outputdir+"/"+species+"_final_annot.gff. Skipping ...")
            file_dict[species].append(outputdir+"/"+species+"_final_annot.gff")

    return file_dict


def filter_te_table(file_dict, outputdir, min_length, max_length):
    for assembly in file_dict.keys():
        nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][7]))
        joined_annot = file_dict[assembly][7]
        if not os.path.exists(outputdir + '/'+''.join(nombre_archivo)+'.bed'):
            df = pd.read_csv(joined_annot, sep="\t", header=None)
            _len = df[4] - df[3]
            keep = _len >= min_length
            if max_length is not None:
                keep &= _len <= max_length
            df_filtered = df[keep]
            df_filtered = df_filtered.reset_index()
            df_filtered['ID'] = df_filtered.apply(
                lambda row: f"{row.name + 1}_{row[8]}_{row[3]}_{row[4]}".replace("Name=", "")
                                                                        .replace("Class=", "-")
                                                                        .replace(";", "-"),
                axis=1
            )

            df_filtered_out = df_filtered[[0, 3, 4, 'ID']]
            df_filtered_out.to_csv(outputdir + '/'+''.join(nombre_archivo)+'.bed', sep="\t", header=False, index=False)
        file_dict[assembly].append(outputdir + '/'+''.join(nombre_archivo)+'.bed')

    return file_dict


def generate_flanking_regions(file_dict, flank_size, outputdir):
    for assembly in file_dict.keys():
        nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][7]))
        bed_path = file_dict[assembly][8]
        if not os.path.exists(f"{outputdir}/Table_{nombre_archivo[0]}_anchors{flank_size}.txt"):
            bed = pd.read_csv(bed_path, sep="\t", header=None,
                              names=["chr", "start", "end", "ID"])
            bed['F1_start'] = bed['start'] - flank_size
            bed['F1_start'] = bed['F1_start'].apply(lambda x: max(x, 1))
            bed['F2_end'] = bed['end'] + flank_size

            f1 = bed[['chr', 'F1_start', 'start', 'ID']].copy()
            f1['ID'] += "_F1"
            f2 = bed[['chr', 'end', 'F2_end', 'ID']].copy()
            f2['ID'] += "_F2"
            f1.columns = ["chr", "start", "end", "ID"]
            f2.columns = ["chr", "start", "end", "ID"]

            combined = pd.concat([f1, f2]).sort_values("ID")

            gff_df = pd.DataFrame()
            gff_df["seqid"] = combined["chr"]
            gff_df["source"] = "Mmd"
            gff_df["type"] = "repeat_region"
            gff_df["start"] = combined["start"]
            gff_df["end"] = combined["end"]
            gff_df["score"] = "."
            gff_df["strand"] = "+"
            gff_df["phase"] = "."
            gff_df["attributes"] = "ID=" + combined["ID"]

            with open(f"{outputdir}/Table_{nombre_archivo[0]}_anchors{flank_size}.txt", 'w') as f:
                f.writelines("##gff-version 3\n#!gff-spec-version 1.21\n#!processor NCBI annotwriter\n")
                gff_df.to_csv(f, sep="\t", index=False, header=False)
        else:
            print(
                "WARNING: Flanked sequence table is already present in: " + outputdir + "/Table_" + nombre_archivo[0] + "_anchors"+str(flank_size)+".txt . Skipping ...")
        file_dict[assembly].append(f"{outputdir}/Table_{nombre_archivo[0]}_anchors{flank_size}.txt")

    return file_dict


def run_liftoff(gff_path, ref_fasta, query_fasta, features_file, output_prefix, output_dir, threads, verbose):
    # Give each liftoff a UNIQUE intermediate dir (derived from the per-job output_prefix) so
    # concurrent jobs (LIFTOFF_MAX_PARALLEL) never collide on the shared reference index / DB that
    # liftoff writes into -dir. The -o outputs are already unique per job.
    inter_dir = f"{output_prefix}_liftoff_tmp"
    os.makedirs(inter_dir, exist_ok=True)
    # Liftoff/minimap2/pysam write sidecar index files (.fai, .mmi) NEXT TO the input fasta paths.
    # When several concurrent jobs share the same underlying reference/query fasta they race on those
    # sidecars -- which silently corrupts one job's index and makes it lift ~nothing (empty polished
    # GFF), not just crash. Point each job at PRIVATE per-job symlinks inside its own inter_dir so
    # every sidecar index is unique to the job. (Symlinks avoid copying the fasta data.)
    tgt_link = os.path.join(inter_dir, "target.fa")
    ref_link = os.path.join(inter_dir, "reference.fa")
    for link, src in ((tgt_link, ref_fasta), (ref_link, query_fasta)):
        if os.path.islink(link) or os.path.exists(link):
            os.remove(link)
        os.symlink(os.path.abspath(src), link)
    command = [
        "liftoff",
        "-g", gff_path,
        "-o", f"{output_prefix}_liffout.gff",
        "-exclude_partial",
        "-dir", inter_dir,
        "-polish",
        "-copies",
        "-f", features_file,
        "-overlap", "1",
        "-p", str(threads),
        tgt_link,
        ref_link
    ]
    output = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if verbose:
        with open(inter_dir + "/Liftoff_log.txt", "w") as logfile:
            logfile.write(output.stdout)
            logfile.write("\n\nSTDERR:\n")
            logfile.write(output.stderr)

    # Surface failures instead of silently producing nothing (the v1 behaviour).
    if output.returncode != 0:
        print(f"WARNING: liftoff failed (return code {output.returncode}) for -g {gff_path}. "
              f"stderr tail: {output.stderr.strip().splitlines()[-1] if output.stderr.strip() else 'n/a'}")

    return f"{output_prefix}_liffout.gff"


def filter_te_insertions(input_path, output_path):
    df = pd.read_csv(input_path, sep="\t", header=None)
    df.columns = ["Chrm", "F1_StartR", "F1_EndR", "F2_StartR", "F2_EndR",
                  "TE_Family", "LengthT", "LengthR", "F1_TE_F2", "State"]

    filtered = df[
        (df["State"] == "Absent") |
        ((df["State"] == "Present") &
         (df["LengthT"] >= df["LengthR"] - 50) &
         (df["LengthT"] <= df["LengthR"] + 50))
    ]
    filtered.to_csv(output_path, sep="\t", index=False, header=False)


def parse_paf_unique_assignments(file_dict, reference, outputdir):
    """
    Procesa el archivo PAF y asigna un unico scaffold a cada cromosoma,
    basado en la mayor longitud de alineamiento (coverage).
    """
    reference_assembly = file_dict[reference][0]
    for assembly in file_dict.keys():
        paf_file = file_dict[assembly][10]
        fasta_in = file_dict[assembly][0]
        anchorfile = file_dict[assembly][9]
        nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][9]))[0]
        if assembly != reference:
            colnames = [
                "query", "qlen", "qstart", "qend", "strand",
                "target", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"
            ]
            df = pd.read_csv(paf_file, sep="\t", header=None, usecols=range(12), names=colnames, comment="#")

            # Calcular longitud alineada
            df["aligned_len"] = df["qend"] - df["qstart"]
            df = df[df["aligned_len"] > 0]

            # Agrupar por pares (target, query) y sumar alineamientos
            pair_cov = (
                df.groupby(["target", "query"])["aligned_len"]
                .sum()
                .reset_index()
            )

            # Para cada cromosoma (target), escoger el scaffold (query) con mayor cobertura
            scaffold_map = (
                pair_cov.sort_values("aligned_len", ascending=False)
                .drop_duplicates("target")
                .set_index("query")["target"]
                .to_dict()
            )

            if not os.path.exists(outputdir+"/"+nombre_archivo+"_chrRenamed.fasta"):
                """
                Cambia los nombres de los scaffolds solo si aparecen en scaffold_map.
                Los demas se mantienen igual.
                """
                renamed = []
                for record in SeqIO.parse(fasta_in, "fasta"):
                    old_id = record.id
                    if old_id in scaffold_map:
                        new_id = scaffold_map[old_id]
                    else:
                        new_id = old_id
                    record.id = new_id
                    record.name = new_id
                    record.description = ""
                    renamed.append(record)

                SeqIO.write(renamed, outputdir+"/"+nombre_archivo+"_chrRenamed.fasta", "fasta")
            else:
                print("WARNING: renamed file already found in " + outputdir+"/"+nombre_archivo+"_chrRenamed.fasta. Skipping....")

            file_dict[assembly].append(outputdir+"/"+nombre_archivo+"_chrRenamed.fasta")

            if not os.path.exists(anchorfile+"_chrRenamed.gff"):
                # To replace Chr names in the bed files
                colnames = ["seqid", "species", "type", "start", "end", "dot", "strand", "dot2", "ID"]
                gff_df = pd.read_csv(anchorfile, sep='\t', comment="#", header=None, names=colnames)

                # Reemplazar los valores de la columna 0 (seqid) segÃºn el diccionario
                gff_df["seqid"] = gff_df["seqid"].map(scaffold_map).fillna(gff_df["seqid"])  # Si no estÃ¡ en el diccionario, conserva el original

                # Guardar el GFF modificado
                gff_df.to_csv(anchorfile+"_chrRenamed.gff", sep="\t", header=False, index=False)
            else:
                print("WARNING: renamed anchor file file already found in " + anchorfile+"_chrRenamed.gff. Skipping....")

            file_dict[assembly].append(anchorfile + "_chrRenamed.gff")

        else:
            file_dict[assembly].append("NoFileNeeded")
            file_dict[assembly].append("NoFileNeeded")

    return file_dict


def run_minimap(file_dict, reference, outputdir, cores):
    reference_assembly = file_dict[reference][0]
    for assembly in file_dict.keys():
        assembly_fasta = file_dict[assembly][0]
        nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][5]))[0]
        if assembly != reference:
            if not os.path.exists(outputdir + '/'+''.join(nombre_archivo)+'.paf'):
                output = subprocess.run(
                    ["minimap2", "-x", "asm5", "-t", str(cores), reference_assembly, assembly_fasta],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                with open(outputdir + '/'+''.join(nombre_archivo)+'.paf', "w") as map_file:
                    map_file.write(output.stdout)
            else:
                print("WARNING: PAF file already found in the path: " + outputdir + '/'+''.join(nombre_archivo)+'.paf. Skipping....')
            file_dict[assembly].append(outputdir + '/'+''.join(nombre_archivo)+'.paf')
        else:
            file_dict[assembly].append("NoFileNeeded")
    return file_dict


# =============================================================================
# TE-transfer detection (v2)
#
# This block replaces the v1 functions extract_te_anchor_pairs /
# build_te_insertion_table / filter_insertions_by_length and the v1
# merge/groupby/collapse helpers. The logic is a faithful port of the
# validated run_pipeline.sh and reproduces its output byte-for-byte. The only
# substitution is the R/dplyr length filter -> an equivalent awk one-liner, so
# BurriTE no longer needs R/dplyr (not part of its conda env).
#
# Each stage is the exact shell command that was validated; it is run via
# subprocess with positional arguments so no shell-escaping of paths is needed.
# Requires only: bash, coreutils (sort/sed/grep/cut/cat), awk and bedtools.
# =============================================================================

# Per species: <gff_polished> <species_label> <out_dir>
#   -> <sp>_F1F2_coordinates.txt  <sp>_TEinsertions.txt
#      Table_<sp>_500filtered.txt  <sp>.bed
_PER_SPECIES_SH = r'''
set -euo pipefail
gff="$1"; sp="$2"; dir="$3"
cd "$dir"

# Steps 1-5: keep base IDs that mapped as exactly one F1_0 + one F2_0
grep -v "#" "$gff" | awk '{print$NF}' | awk -F";" '{print$NF}' | grep -E "F1_0|F2_0" \
  | sed 's/copy_num_//g' | sed 's/_F1_/\t/g' | sed 's/_F2_/\t/g' | awk '{print$1}' \
  | sort | uniq -c | awk '$1==2{print $2}' > "_basekeys_${sp}.txt"
awk '{print $0"_F1_0"; print $0"_F2_0"}' "_basekeys_${sp}.txt" | sort -n > "_tempF1F2_${sp}.txt"

# Step 7: per-flank coordinates keyed by copy_num_ID
grep -v "#" "$gff" | sed 's/;copy_num_ID/\tID/g' | awk '{print$10"\t"$1"\t"$4"\t"$5"\t"$10}' > "_tabtemp_${sp}.txt"

# Step 9: attach coordinates to each F1_0 / F2_0 (NA if unresolved)
awk 'NR==FNR{a[$1]=$0;next}; {print $1, $1 in a?a[$1]:"NA"}' "_tabtemp_${sp}.txt" "_tempF1F2_${sp}.txt" \
  | cut -f2- | sed 's/_ID=/_/g' > "${sp}_F1F2_coordinates.txt"
rm -f "_tempF1F2_${sp}.txt" "_tabtemp_${sp}.txt" "_basekeys_${sp}.txt"

# Pair F1<->F2 by shared ID; build short name + length; call Present/Absent
awk -F"\t" -v sp="$sp" '
  NF==4 && $4 ~ /_F[12]_0$/ {
    base=$4; sub(/_F[12]_0$/,"",base);
    if ($4 ~ /_F1_0$/) { f1c[base]=$1; f1s[base]=$2; f1e[base]=$3; seenF1[base]=1 }
    else               { f2s[base]=$2; f2e[base]=$3; seenF2[base]=1 }
    if(!(base in order)){ order[base]=++n; keys[n]=base }
  }
  END{
    for(i=1;i<=n;i++){ b=keys[i];
      if(!(seenF1[b] && seenF2[b])) continue;
      m=split(b,a,"_"); len=a[m]-a[m-1];
      name=b; sub(/_[0-9]+_[0-9]+$/,"",name);
      fam="";
      if (name ~ /--/ && match(name, /^ID=[0-9]+/)) {
        pre=substr(name,1,RLENGTH);
        rest=substr(name,RLENGTH+1);        # "_TE<rmnum>-<family>--<class>"
        sf=rest; sub(/^.*--/,"",sf);          # class (after the "--" Name/Class boundary)
        fam=rest; sub(/--.*$/,"",fam); sub(/^_TE[0-9]+-/,"",fam);   # family (Name, before "--")
        name=pre"_"sp"_"sf;
      }
      gsub(/\//,"_",name); gsub(/\//,"_",fam);
      if (fam!="") name=name"#"fam;           # carry family as a #-delimited suffix
      gap=f2s[b]-f1e[b];
      state=(gap>=0 && gap<=10)?"Absent":"Present";
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s-%s \t%s\n", f1c[b],f1s[b],f1e[b],f2s[b],f2e[b],name,len,gap,f1e[b],f2s[b],state;
    }
  }
' "${sp}_F1F2_coordinates.txt" > "${sp}_TEinsertions.txt"

# Length filter (R/dplyr replaced by equivalent awk): keep Absent or |LengthT-LengthR|<=20% of LengthT
awk -F"\t" 'BEGIN{OFS="\t"; print "Chrm","F1_StartR","F1_EndR","F2_StartR","F2_EndR","TE_Family","LengthT","LengthR","F1_TE_F2","State"}
  $10=="Absent" || ($10=="Present" && $8>=$7-0.20*$7 && $8<=$7+0.20*$7)' "${sp}_TEinsertions.txt" > "Table_${sp}_500filtered.txt"

# Per-species BED: Chrm  F1_End  F2_Start  ID|LengthT
#   The |value carries LengthT (column 7 = the real TE length), NOT LengthR (the mapped
#   flank gap). For Present calls LengthT ~ LengthR, but for Absent calls the gap is ~0-10
#   while the TE is large; carrying LengthT keeps absence-derived loci correctly sized so
#   the collapse and step-14 genotyping use the true TE length instead of the tiny gap.
awk '{print$1"\t"$3"\t"$4"\t"$6"|"$7}' "Table_${sp}_500filtered.txt" | sed '1d' > "${sp}.bed"
'''

# Merge all per-species beds: <bed_dir> <out_bed>
_MERGE_SH = r'''
set -euo pipefail
cat "$1"/*.bed | sort -k1,1 -k2,2n | awk -F"\t" 'BEGIN{OFS="\t"}{
  split($4, a, "|"); len=a[2];
  split(a[1], h, "#"); core=h[1];          # drop the optional #family suffix (groupby = superfamily only)
  n=split(core, b, "_");
  species=b[2];
  superfam=b[3]; for(i=4;i<=n;i++) superfam=superfam"_"b[i];
  print $1,$2,$3,superfam,len,species
}' > "$2"
'''

# bedtools groupby on superfamily + rebuild one interval per group: <merged> <out>
_GROUPBY_SH = r'''
set -euo pipefail
bedtools groupby -i "$1" -g 4 -c 1,2,3,4,5,6 -o collapse,collapse,collapse,distinct,collapse,distinct | awk 'BEGIN { FS=OFS="\t" } {split($2, col2, ","); split($3, col3, ","); split($4, col4, ","); unique_col2 = col2[1]; min_col3 = col3[1]; for (i = 2; i <= length(col3); i++) {if (col3[i] < min_col3) { min_col3 = col3[i];}} max_col4 = col4[1]; for (i = 2; i <= length(col4); i++) {if (col4[i] > max_col4) {max_col4 = col4[i];}} print unique_col2"\t" min_col3"\t" max_col4"\t"$5"\t" $6"\t"$7;}' > "$2"
'''

# Length-consistency filter + final collapsed BED with FM<index> name: <grouped> <out>
# Input groupby.bed cols: Chrm start end superfamily length(list=LengthT) species
#   1) keep groups whose member LengthT agree within 50 bp (real-size consistency)
#   2) label each surviving group FM<idx>_<superfamily>|<mean LengthT>  (the real TE size)
# The old coord-span filter (drop if span <= 0.5*mean) is removed: now that the length is
# the real TE size (not the flank gap), absence-derived loci legitimately have a tiny
# coord span with a large length, and must be kept rather than dropped.
_COLLAPSE_SH = r'''
set -euo pipefail
awk -F"\t" '{split($5, v, ","); min = max = v[1]; for (i in v) { if (v[i] < min) min = v[i]; if (v[i] > max) max = v[i]; } print $0 "\t" ((max - min > 50) ? "FALSE" : "TRUE"); }' "$1" \
  | grep "TRUE" \
  | awk -F"\t" '{split($5, v, ","); sum = 0; for (i in v) sum += v[i]; mean = int(sum / length(v) + 0.5); print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"mean; }' \
  | cat -n \
  | awk '{print $2"\t"$3"\t"$4"\tFM"$1"_"$5"|"$7}' > "$2"
'''


def _run_sh(script, *args):
    """Run an embedded shell block with positional args ($1, $2, ...)."""
    subprocess.run(["bash", "-c", script, "bash", *map(str, args)], check=True)


def te_transfer_per_species(gff_polished, species_label, out_dir):
    """Port of the per-species run_pipeline.sh block (validated byte-for-byte).

    Reads a Liftoff *_polished GFF whose TEs are carried as F1_0/F2_0 flank
    markers in copy_num_ID, and writes into out_dir:
        <sp>_F1F2_coordinates.txt   per-flank coords
        <sp>_TEinsertions.txt       paired table with Present/Absent calls
        Table_<sp>_500filtered.txt  length-filtered final table
        <sp>.bed                    Chrm F1_End F2_Start ID|LengthR
    Returns the path to <sp>.bed.
    """
    # The shell block does `cd "$dir"`, so the input GFF must be an absolute path
    # (otherwise a relative path stops resolving after the cd). out_dir is made
    # absolute too so the returned bed path is stable regardless of CWD.
    gff_polished = os.path.abspath(gff_polished)
    out_dir = os.path.abspath(out_dir)
    _run_sh(_PER_SPECIES_SH, gff_polished, species_label, out_dir)
    return os.path.join(out_dir, f"{species_label}.bed")


def transfer_anchors(chromosome_list, file_dict, reference_name, cores, verbose, outputdir):
    with open(chromosome_list) as chrm_file:
        chromosomes = [line.strip() for line in chrm_file]

    # Create a folder to save the final bed files
    create_output_folders(outputdir+"/bed_files")

    # --- Phase 1: liftoff the anchors per (chromosome, assembly) -> polished GFFs ---
    # Prepare per-(chrm, assembly) inputs serially (fast I/O), collect the liftoff jobs, then run
    # them concurrently in a pool of LIFTOFF_MAX_PARALLEL workers (each liftoff uses its own -dir).
    polished_by_assembly = defaultdict(list)
    liftoff_jobs = []
    per_job_threads = max(1, int(cores) // LIFTOFF_MAX_PARALLEL)
    for chrm in chromosomes:
        print(f"INFO: Processing chromosome: {chrm}")
        # To create a separate fasta file for each desired chromosome in the reference
        assembly_fasta = file_dict[reference_name][0]
        seqiter = SeqIO.parse(assembly_fasta, "fasta")
        nombre_archivo = os.path.splitext(os.path.basename(file_dict[reference_name][0]))[0]
        with open(outputdir + "/" + nombre_archivo + "_" + chrm + ".fasta", 'w') as output_file:
            SeqIO.write((seq for seq in seqiter if seq.id == chrm), output_file, "fasta")

        for assembly in file_dict.keys():
            if assembly == reference_name:
                continue
            polished_gff = f"{outputdir}/{assembly}_{chrm}_liffout.gff_polished"
            if not os.path.exists(polished_gff):
                # To create a separate fasta file for each desired chromosome in the alternative assemblies
                assembly_fasta = file_dict[assembly][11]
                seqiter = SeqIO.parse(assembly_fasta, "fasta")
                nombre_archivo = os.path.splitext(os.path.basename(file_dict[assembly][0]))[0]
                with open(outputdir + "/" + nombre_archivo + "_" + chrm + ".fasta", 'w') as output_file:
                    num_seqs = SeqIO.write((seq for seq in seqiter if seq.id == chrm), output_file, "fasta")
                    if num_seqs == 0:
                        print(f"ERROR: No sequences found named {chrm} in file {assembly_fasta}")
                        sys.exit(0)

                query_fasta = outputdir+ "/"+nombre_archivo+"_"+chrm+".fasta"
                nombre_archivo_ref = os.path.splitext(os.path.basename(file_dict[reference_name][0]))[0]
                ref_fasta = outputdir+ "/"+nombre_archivo_ref+"_"+chrm+".fasta"

                chrm_gff = f"{outputdir}/{assembly}_anchors{flank_size}_{chrm}.gff"
                colnames = ["seqid", "species", "type", "start", "end", "dot", "strand", "dot2", "ID"]
                gff_df = pd.read_csv(file_dict[assembly][12], sep='\t', comment="#", header=None, names=colnames)
                gff_df_chr = gff_df[gff_df["seqid"] == chrm]
                gff_df_chr.to_csv(chrm_gff, sep="\t", index=False, header=False)

                features_file = outputdir + "/../04_graffite_liftover/features_file.txt"
                # Pre-index the (shared) reference and query fastas serially so concurrent jobs don't
                # race building the .fai (see ensure_faidx).
                ensure_faidx(ref_fasta)
                ensure_faidx(query_fasta)
                # run_liftoff is called with -polish, which writes "<o>_liffout.gff_polished"
                liftoff_jobs.append((chrm_gff, ref_fasta, query_fasta, features_file,
                                     f"{outputdir}/{assembly}_{chrm}", outputdir, per_job_threads, verbose))
            else:
                print(f"WARNING: Polished liftoff already present: {polished_gff}. Skipping liftoff ...")
            polished_by_assembly[assembly].append(polished_gff)

    if liftoff_jobs:
        print(f"INFO: running {len(liftoff_jobs)} step-10 liftoff jobs, up to {LIFTOFF_MAX_PARALLEL} "
              f"in parallel ({per_job_threads} internal threads each)")
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(LIFTOFF_MAX_PARALLEL, len(liftoff_jobs))) as ex:
            list(ex.map(lambda a: run_liftoff(*a), liftoff_jobs))

    # --- Phase 2: per-species TE-transfer detection on the polished GFFs ---
    #     (validated port of run_pipeline.sh; uses the *polished* GFF, all chromosomes combined)
    for assembly, polished_list in polished_by_assembly.items():
        species_bed = f"{outputdir}/bed_files/{assembly}.bed"
        if os.path.exists(species_bed):
            print(f"WARNING: Filtered bed file already present in : {species_bed}. Skipping ...")
            continue
        # The TE family token embedded in the bed is ID=<num>_<species>_<superfamily>,
        # which the merge step (step 11) parses by splitting on "_". The species label
        # must therefore be underscore-free, otherwise an assembly name like "C_arabica"
        # leaks into the superfamily field (-> "arabica_LTR") and corrupts the grouping
        # key so cross-species loci no longer collapse. Stripping "_" reproduces the
        # validated convention (C_arabica -> Carabica).
        species_tag = assembly.replace("_", "")
        combined_gff = f"{outputdir}/{assembly}_combined_liffout.gff_polished"
        with open(combined_gff, "w") as out:
            for p in polished_list:
                with open(p) as fh:
                    shutil.copyfileobj(fh, out)
        bed_path = te_transfer_per_species(combined_gff, species_tag, outputdir)
        shutil.copy(bed_path, species_bed)

    return outputdir+"/bed_files/"


def merge_beds_to_manaus_nhv(input_bed_dir, output_bed):
    """Merge all per-species beds, sort by position and split the ID column into
    Chrm start end superfamily length species (validated run_pipeline.sh logic).

    NOTE: replaces the v1 regex parser, which silently dropped rows whose family
    used the '--' separator (it only matched a literal '_LTR'/'_DNA'/...).
    """
    _run_sh(_MERGE_SH, input_bed_dir, output_bed)


def run_bedtools_groupby(input_bed, output_bed):
    """Group rows by superfamily (col 4) with bedtools and rebuild one interval
    per group: Chrm(first) min(start) max(end) superfamily length(list) species.

    NOTE: the v2 grouping key is the superfamily (built in merge_beds_to_manaus_nhv),
    so groupby actually clusters the same insertion seen across species; the v1
    key was unique per row, making groupby a no-op.
    """
    _run_sh(_GROUPBY_SH, input_bed, output_bed)


def collapse_and_filter_grouped_bed(grouped_bed, output_bed):
    """Keep length-consistent groups (member lengths within 50 bp), drop tiny
    spans, and write the final collapsed BED with an FM<index> name
    (validated run_pipeline.sh logic)."""
    _run_sh(_COLLAPSE_SH, grouped_bed, output_bed)


def _median_int(vals):
    """Integer median matching the awk median() in build_incremental_collapse.sh."""
    b = sorted(vals); n = len(b)
    if n % 2:
        return b[(n - 1) // 2]
    return int((b[n // 2 - 1] + b[n // 2]) / 2)


def incremental_collapse(bed_paths, samples, out_dir, prefix,
                         lenthr=100, minlen=0, bedmin=100, merge_by="superfamily"):
    """Incremental locus+length TE collapse (alternative to the superfamily groupby).

    merge_by (--merge_by): the clustering key. "superfamily" (default) groups the Class-level
    token (e.g. DNA_RC); "family" groups the finer Name-level token (e.g. ATREP13), which is
    stricter (splits a superfamily locus into its distinct families) and better suited to
    within-species comparisons. Each cluster records BOTH: its superfamily and the SET of
    member families (so a superfamily locus that merged several families is visible). Loci are
    given opaque ids LOC<n>; family/superfamily are carried as separate columns, not baked
    into the id.

    Python port of build_incremental_collapse.sh. A per-individual interval joins an
    existing consensus cluster ONLY if all hold: (1) same chromosome, (2) same TE
    family, (3) coordinate overlap, (4) the cluster's internal length range stays
    <= lenthr. Otherwise it seeds a new cluster, so length variants and private loci
    are kept (unlike the superfamily groupby, which merges them away). Each cluster is
    labelled shared (>=2 individuals) or singleton (==1).

    The .tsv and dictionary carry a Nested column (TRUE if the consensus TE is fully
    contained within another consensus TE of a different grouping key = true TE-in-TE
    nesting; same-key length variants stay FALSE).

    Writes into out_dir:
        <prefix>_collapsed.tsv     full catalogue (... superfamily family n_families families samples Nested)
        <prefix>_shared.bed        clusters seen in >=2 individuals
        <prefix>_singletons.bed    clusters seen in 1 individual
        <prefix>_collapsed.bed     4-col consensus BED, median length >= bedmin (downstream input)
        <prefix>_TE_collapsed.txt  the id list from <prefix>_collapsed.bed
        <prefix>_dictionary.tsv    consensus TE -> original per-individual TE ids (+ Nested)
    Returns the path to <prefix>_collapsed.bed.

    Family is parsed as in merge_beds_to_manaus_nhv: from ID=<num>_<sample>_<family>,
    family = everything after the 2nd underscore (so it is independent of the file name).
    """
    # 1) tagged intervals: (chrom, start, end, family, sample, original_id)
    rows = []
    for path, sample in zip(bed_paths, samples):
        if not os.path.exists(path):
            print(f"WARNING: incremental collapse: missing {path}; skipping {sample}.")
            continue
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                chrom, s, e = parts[0], int(parts[1]), int(parts[2])
                if e < s:
                    s, e = e, s
                oid = parts[3]              # ID=<num>_<sample>_<superfam>[#<family>]|<LengthT>
                name, _, lp = oid.partition("|")
                # length is the real TE size (LengthT) from the |value, not the coord width
                reclen = int(lp) if lp.lstrip("-").isdigit() else (e - s)
                if reclen < minlen:
                    continue
                core, _, fam = name.partition("#")           # superfam part / family part
                flds = core.split("_")
                superfam = "_".join(flds[2:]) if len(flds) > 2 else core
                family = fam if fam else superfam            # fall back if no #family present
                key = family if merge_by == "family" else superfam
                rows.append((chrom, s, e, key, superfam, family, sample, oid, reclen))
    # sort matches `sort -k1,1 -k<key>,<key> -k2,2n` + GNU last-resort whole-line tiebreak
    rows.sort(key=lambda r: (r[0], r[3], r[1],
                             f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[6]}\t{r[7]}"))

    # 2) cluster (same chrom + same key + overlap + internal length range <= lenthr)
    clusters, cur = [], None
    for chrom, s, e, key, superfam, family, sample, oid, reclen in rows:
        tl = reclen
        if cur and chrom == cur["chrom"] and key == cur["key"] and s <= cur["end"] \
                and (max(tl, cur["lmax"]) - min(tl, cur["lmin"])) <= lenthr:
            cur["end"] = max(cur["end"], e)
            cur["lmin"] = min(cur["lmin"], tl); cur["lmax"] = max(cur["lmax"], tl)
            cur["lengths"].append(tl); cur["samples"].add(sample)
            cur["families"].add(family)
            if oid not in cur["oids_set"]:
                cur["oids_set"].add(oid); cur["oids"].append(oid)
        else:
            if cur:
                clusters.append(cur)
            cur = {"chrom": chrom, "start": s, "end": e, "key": key, "superfam": superfam,
                   "families": {family}, "lengths": [tl], "lmin": tl, "lmax": tl,
                   "samples": {sample}, "oids": [oid], "oids_set": {oid}}
    if cur:
        clusters.append(cur)

    # 2b) NESTED flag: TRUE if a consensus TE is fully contained within another consensus TE
    #     of a DIFFERENT grouping key (true TE-in-TE nesting across families); same-key length
    #     variants stay FALSE. Sweep sorted by chrom, start ASC, end DESC: a host processed
    #     before the current locus has start <= current.start, so if any OTHER key's running
    #     max end >= current.end, the current locus is nested within it. (Port of v2's step 2b.)
    cur_chrom, maxend = None, {}
    for cl in sorted(clusters, key=lambda c: (c["chrom"], c["start"], -c["end"])):
        if cl["chrom"] != cur_chrom:
            cur_chrom, maxend = cl["chrom"], {}
        e, k, nz = cl["end"], cl["key"], False
        for g, ge in maxend.items():
            if g != k and ge >= e:
                nz = True
                break
        cl["nested"] = nz
        if e > maxend.get(k, -1):
            maxend[k] = e

    # 3) write outputs
    p_tsv  = os.path.join(out_dir, f"{prefix}_collapsed.tsv")
    p_shar = os.path.join(out_dir, f"{prefix}_shared.bed")
    p_sing = os.path.join(out_dir, f"{prefix}_singletons.bed")
    p_bed  = os.path.join(out_dir, f"{prefix}_collapsed.bed")
    p_ids  = os.path.join(out_dir, f"{prefix}_TE_collapsed.txt")
    p_dict = os.path.join(out_dir, f"{prefix}_dictionary.tsv")

    bed_rows = []
    with open(p_tsv, "w") as tsv, open(p_shar, "w") as shar, \
            open(p_sing, "w") as sing, open(p_dict, "w") as dic:
        tsv.write("Chrm\tstart\tend\tLOC_ID\tn_individuals\tclass\tlen_min\tlen_max\t"
                  "superfamily\tfamily\tn_families\tfamilies\tsamples\tNested\n")
        dic.write("Chrm\tstart\tend\tLOC_ID\tn_individuals\toriginal_ids\tNested\n")
        for idx, cl in enumerate(clusters, 1):
            n = len(cl["samples"]); cls = "singleton" if n == 1 else "shared"
            med = _median_int(cl["lengths"])
            loc = f"LOC{idx}|{med}"                          # opaque locus id (+ median length)
            fams = sorted(cl["families"]); nf = len(fams)
            fam_field = fams[0] if nf == 1 else "mixed"      # single family, else "mixed"
            nested = "TRUE" if cl["nested"] else "FALSE"
            row4 = f"{cl['chrom']}\t{cl['start']}\t{cl['end']}\t{loc}"
            tsv.write(f"{row4}\t{n}\t{cls}\t{cl['lmin']}\t{cl['lmax']}\t{cl['superfam']}\t"
                      f"{fam_field}\t{nf}\t{';'.join(fams)}\t{','.join(cl['samples'])}\t{nested}\n")
            (shar if cls == "shared" else sing).write(row4 + "\n")
            dic.write(f"{cl['chrom']}\t{cl['start']}\t{cl['end']}\t{loc}\t{n}\t{','.join(cl['oids'])}\t{nested}\n")
            if med >= bedmin:
                bed_rows.append((cl["chrom"], cl["start"], cl["end"], loc))

    bed_rows.sort(key=lambda r: (r[0], r[1]))            # position-sorted 4-col BED
    with open(p_bed, "w") as bed, open(p_ids, "w") as ids:
        for chrom, s, e, loc in bed_rows:
            bed.write(f"{chrom}\t{s}\t{e}\t{loc}\n")
            ids.write(loc + "\n")

    shared_n = sum(1 for cl in clusters if len(cl["samples"]) >= 2)
    print(f"INFO: incremental collapse -> {len(clusters)} consensus TEs "
          f"({shared_n} shared, {len(clusters)-shared_n} singleton); "
          f"{len(bed_rows)} kept in {os.path.basename(p_bed)} (median >= {bedmin} bp).")
    return p_bed


def build_reference_bed(file_dict, reference_name, min_length, max_length, out_bed_dir):
    """Emit the reference genome's OWN TE annotation as one more per-genome BED so
    reference(-private) loci enter the step-11 consensus catalogue (--include_reference).

    The reference TEs are already in reference coordinates (they need no ALT->reference
    liftoff, unlike the query loci), so this simply reformats the reference's joined
    RepeatMasker+GraffiTE annotation (slot 7, <ref>_final_annot.gff) into the exact BED
    that merge_beds_to_manaus_nhv / incremental_collapse consume:
        Chrm  TE_start  TE_end  ID=<n>_<reftag>_<superfamily>|<LengthT>
    where reftag = reference_name with '_' stripped (kept underscore-free so the merge's
    split-on-'_' puts it in the species field, mirroring the species_tag convention), the
    superfamily is the Class= attribute with '/'->'_', and LengthT = end - start (the real
    reference TE size, the analogue of OneCode column 7). The same min/max length filter as
    filter_te_table is applied. Returns the path to <reftag>.bed.
    """
    ref_tag = reference_name.replace("_", "")
    out_bed = os.path.join(out_bed_dir, f"{ref_tag}.bed")
    joined_annot = file_dict[reference_name][7]
    df = pd.read_csv(joined_annot, sep="\t", header=None, comment="#",
                     names=["seqid", "source", "type", "start", "end",
                            "score", "strand", "phase", "attributes"])
    n = 0
    with open(out_bed, "w") as out:
        for _, row in df.iterrows():
            start, end = int(row["start"]), int(row["end"])
            if end < start:
                start, end = end, start
            length = end - start
            if length < min_length or (max_length is not None and length > max_length):
                continue
            m = re.search(r"Class=([^;\s]+)", str(row["attributes"]))
            superfam = (m.group(1) if m else "Unknown").replace("/", "_")
            # family = the Name= attribute (the token before '#' in the library), carried as a
            # #-delimited suffix so --merge_by family can group the reference too.
            mf = re.search(r"Name=([^;\s]+)", str(row["attributes"]))
            family = (mf.group(1) if mf else "").replace("/", "_")
            n += 1
            tag = f"ID={n}_{ref_tag}_{superfam}" + (f"#{family}" if family else "")
            out.write(f"{row['seqid']}\t{start}\t{end}\t{tag}|{length}\n")
    print(f"INFO: --include_reference: wrote {n} reference TE loci to {out_bed}")
    return out_bed


def evaluate_reference_transfer(collapsed_bed, output_path, flank_size):
    """Genotype the reference genome itself at every consensus FM locus (--include_reference).

    The consensus flanks are in reference coordinates, so lifting them back onto the
    reference is the identity and the reference F1_End/F2_Start are just the FM locus
    start/end; the reference gap is therefore end - start. The Present/Absent/Ambiguous
    definition is identical to evaluate_te_transfer, but there is no NA (both flanks are by
    construction present in the reference). Writes the same 10-column <ref>_Transfers_TE.txt
    that build_presence_absence_matrix reads, so the reference becomes an ordinary column.
    """
    with open(collapsed_bed) as f, open(output_path, "w") as out:
        out.write("Chrm\tF1_StartR\tF1_EndR\tF2_StartR\tF2_EndR\tID\tTE_Length\tGap\tState\tUnmappable_flanks\n")
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            chrm, start, end, tag = line.split("\t")[:4]
            start, end = int(start), int(end)
            name, _, span = tag.partition("|")
            L = int(span) if span.lstrip("-").isdigit() else None
            gap = end - start
            if L is None:
                state = "NA"
            elif 0 <= gap <= 10:
                state = "Absent"
            elif abs(L - gap) <= 0.20 * L:      # Present tolerance scales with TE length (20% of L)
                state = "Present"
            else:
                state = "Ambiguous"
            Lout = L if L is not None else "NA"
            out.write(f"{chrm}\t{start - flank_size}\t{start}\t{end}\t{end + flank_size}\t"
                      f"{name}\t{Lout}\t{gap}\t{state}\t-\n")


def generate_flanking_gffs(collapsed_bed, output_dir, flank_size, chromosome_list):
    # v2 collapsed BED has 4 columns: Chrm  start  end  FM<n>_<superfamily>|<span>
    # (the validated run_pipeline.sh format). The feature label is the FM<n>_...
    # part before the "|"; FM<n> is a running index so labels stay unique.
    df = pd.read_csv(collapsed_bed, sep="\t", header=None,
                     names=["chr", "start", "end", "name"])
    df["label"] = df["name"].str.split("|").str[0]

    header = "##gff-version 3\n"

    os.makedirs(output_dir, exist_ok=True)
    for chrm in chromosome_list:
        df_chr = df[df["chr"] == chrm]
        if df_chr.empty:
            continue

        f1 = df_chr.copy()
        f1["start"] = (f1["start"] - flank_size).clip(lower=1)
        f1["end"] = df_chr["start"]
        f1["ID"] = f1["label"] + "_F1"

        f2 = df_chr.copy()
        f2["start"] = df_chr["end"]
        f2["end"] = df_chr["end"] + flank_size
        f2["ID"] = f2["label"] + "_F2"

        df_out = pd.concat([f1, f2])
        df_out["source"] = "Mmd"
        df_out["type"] = "repeat_element"
        df_out["score"] = "."
        df_out["strand"] = "+"
        df_out["phase"] = "."
        df_out["attributes"] = "ID=" + df_out["ID"]

        df_out = df_out[["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]]
        gff_path = os.path.join(output_dir, f"FM_{chrm}.gff")

        with open(gff_path, 'w') as out:
            out.write(header)
            df_out.to_csv(out, sep="\t", index=False, header=False)

        features_path = os.path.join(output_dir, f"{chrm}_gff_features.txt")
        df_out["type"].to_csv(features_path, index=False, header=False)


def extract_chromosome_fasta(full_fasta, chromosome_list, output_dir, genome_prefix):
    records = SeqIO.to_dict(SeqIO.parse(full_fasta, "fasta"))
    os.makedirs(output_dir, exist_ok=True)
    for chrm in chromosome_list:
        if chrm in records:
            output_fasta = os.path.join(output_dir, f"{genome_prefix}_{chrm}.fasta")
            with open(output_fasta, "w") as out_f:
                SeqIO.write(records[chrm], out_f, "fasta")


def parse_liftoff_polished(gff_path):
    mapping = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            # Key on copy_num_ID, which carries the _F1_0 / _F2_0 primary-copy suffix
            # that evaluate_te_transfer looks up (the bare ID= lacks the _0). Fall back
            # to ID= if copy_num_ID is absent.
            match = re.search(r'copy_num_ID=([^;\s]+)', parts[8]) or re.search(r'ID=([^;]+)', parts[8])
            if match:
                te_id = match.group(1)
                mapping[te_id] = (chrom, start, end)
    return mapping


def evaluate_te_transfer(te_list_path, liftoff_gff_path, output_path):
    """Genotype each consensus FM locus in one assembly from its lifted flanks.

    Uses the same Present/Absent definition as steps 10-11 / run_pipeline.sh, with the
    FM span (from the collapsed BED, e.g. FM1_LTR|723 -> 723) as the TE length:
        gap = F2_start - F1_end
        Absent     if 0 <= gap <= 10        (flanks adjacent -> TE not there)
        Present    if |TE_length - gap| <= 0.20*TE_length  (gap ~ TE length -> TE there)
        Ambiguous  otherwise                 (incl. negative gaps and mis-mapped flanks)
        NA         if a flank didn't lift, or the two flanks landed on different chromosomes

    The Unmappable_flanks column records why a locus is NA: "F1", "F2" or "F1,F2"
    for flank(s) that failed to lift, "diff_chrom" if the flanks landed on different
    chromosomes, or "-" when both flanks lifted to the same chromosome (a real call).
    """
    # FM<n>_<superfamily> -> TE length, read from the collapsed-BED last column FM..|<span>
    te_len, order = {}, []
    with open(te_list_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, _, span = line.split("\t")[-1].partition("|")
            te_len[name] = int(span) if span.lstrip("-").isdigit() else None
            order.append(name)

    anchor_map = parse_liftoff_polished(liftoff_gff_path)

    def calculate_result(te_id):
        L = te_len.get(te_id)
        Lout = L if L is not None else "NA"
        f1 = anchor_map.get(te_id + "_F1_0")
        f2 = anchor_map.get(te_id + "_F2_0")

        if not f1 or not f2:
            # one or both flanks failed to lift; keep whatever coords we have
            unmap = "F1,F2" if (not f1 and not f2) else ("F1" if not f1 else "F2")
            chrm = (f1 or f2 or ("NA",))[0]
            s1, e1 = (f1[1], f1[2]) if f1 else ("NA", "NA")
            s2, e2 = (f2[1], f2[2]) if f2 else ("NA", "NA")
            return f"{chrm}\t{s1}\t{e1}\t{s2}\t{e2}\t{te_id}\t{Lout}\tNA\tNA\t{unmap}\n"

        chrom1, start1, end1 = f1
        chrom2, start2, end2 = f2
        if chrom1 != chrom2:
            return f"{chrom1}\t{start1}\t{end1}\t{start2}\t{end2}\t{te_id}\t{Lout}\tNA\tNA\tdiff_chrom\n"
        gap = start2 - end1
        if L is None:
            state = "NA"
        elif 0 <= gap <= 10:
            state = "Absent"
        elif abs(L - gap) <= 0.20 * L:      # Present tolerance scales with TE length (20% of L)
            state = "Present"
        else:
            state = "Ambiguous"
        return f"{chrom1}\t{start1}\t{end1}\t{start2}\t{end2}\t{te_id}\t{Lout}\t{gap}\t{state}\t-\n"

    with open(output_path, 'w') as out:
        out.write("Chrm\tF1_StartR\tF1_EndR\tF2_StartR\tF2_EndR\tID\tTE_Length\tGap\tState\tUnmappable_flanks\n")
        for te_id in order:
            out.write(calculate_result(te_id))


def build_presence_absence_matrix(collapsed_bed, assemblies, final_annots_dir, output_path,
                                  contrib=None, annot_label="Present:annot", loc_meta=None):
    """Join the per-assembly <assembly>_Transfers_TE.txt genotype tables into one
    presence/absence matrix, one row per consensus locus (reference coords).

    loc_meta (incremental): {LOC<n> -> (family, superfamily)}. When given, the locus id is the
    opaque LOC<n> and family + superfamily are written as their own columns:
        Chrm  start  end  LOC_ID  family  superfamily  TE_length  <assembly1> ...
    When None (groupby legacy), the id/superfamily are parsed from the FM<n>_<superfamily> name:
        Chrm  start  end  FM_ID  superfamily  TE_length  <assembly1>  <assembly2> ...
    Each assembly column holds that assembly's State (Present/Absent/Ambiguous/NA).
    NA cells are annotated with the unmappable flank(s) as "NA:F1", "NA:F2",
    "NA:F1,F2" or "NA:diff_chrom" so missing-data reasons are visible in the matrix.
    A locus absent from an assembly's table (shouldn't happen) defaults to NA.

    Annotation-Present override (always on): contrib is {FM<n>_<superfamily> ->
    set(genome names)} of the genomes that CONTRIBUTED an annotated TE to each consensus
    locus (the step-11 collapse membership, built by build_contributor_map for either merge
    method). A contributing genome has the TE by direct annotation, which is more reliable
    than the reference->genome flank round-trip; so when its liftback genotype is not already
    Present (Absent/Ambiguous/NA, e.g. a short repetitive element whose flanks mis-lift), its
    cell is overridden to `annot_label` ("Present:annot"). Cells already genotyped Present are
    left untouched. Anything startswith "Present" is present. (contrib=None disables it.)
    """
    # per-assembly: FM<n>_<superfamily> -> State (NA annotated with unmappable flanks)
    states = {}
    for sp in assemblies:
        path = os.path.join(final_annots_dir, f"{sp}_Transfers_TE.txt")
        d = {}
        if os.path.exists(path):
            with open(path) as f:
                next(f, None)  # skip header
                for line in f:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) >= 10:
                        state, unmap = cols[8], cols[9]          # State , Unmappable_flanks
                        d[cols[5]] = f"NA:{unmap}" if state == "NA" else state
                    elif len(cols) >= 9:
                        d[cols[5]] = cols[8]
        else:
            print(f"WARNING: matrix build: missing {path}; column filled with NA.")
        states[sp] = d

    overridden = 0            # cells rescued to Present:annot
    override_from = defaultdict(int)   # what the liftback had called before the override
    with open(collapsed_bed) as f, open(output_path, "w") as out:
        if loc_meta is not None:
            out.write("Chrm\tstart\tend\tLOC_ID\tfamily\tsuperfamily\tTE_length\t"
                      + "\t".join(assemblies) + "\n")
        else:
            out.write("Chrm\tstart\tend\tFM_ID\tsuperfamily\tTE_length\t"
                      + "\t".join(assemblies) + "\n")
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            chrm, start, end, tag = line.split("\t")[:4]
            name, _, span = tag.partition("|")        # LOC<n> (incremental) or FM<n>_<superfam>
            if loc_meta is not None:
                family, superfam = loc_meta.get(name, ("NA", "NA"))
                fixed = [chrm, start, end, name, family, superfam, span]
            else:
                fm_id, _, superfam = name.partition("_")   # FM<n> , <superfamily>
                fixed = [chrm, start, end, fm_id, superfam, span]
            cells = []
            members = contrib.get(name, ()) if contrib else ()
            for sp in assemblies:
                st = states[sp].get(name, "NA")
                if sp in members and not str(st).startswith("Present"):
                    override_from[st.split(":")[0]] += 1
                    overridden += 1
                    st = annot_label
                cells.append(st)
            out.write("\t".join(fixed + cells) + "\n")

    if contrib is not None:
        print(f"INFO: annotation-Present override: {overridden} cells set to {annot_label} "
              f"from a contributing genome's direct annotation (was {dict(override_from)}).")


def _gff_attr(attrs, key):
    """Value of key=... in a GFF attributes string, or None."""
    for kv in attrs.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def _fasta_headers(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def _build_seqmap(orig_fasta, renamed_fasta):
    """{original scaffold id -> reference-chromosome id}. The step-9 chr rename is a pure 1:1
    seqid rename (coordinates unchanged; scaffolds with no reference chr keep their own name), so
    pairing the original and renamed FASTA headers by position recovers the map. Empty dict (=>
    identity) when either FASTA is absent (e.g. the reference, already in reference coords)."""
    if not (orig_fasta and renamed_fasta
            and os.path.exists(orig_fasta) and os.path.exists(renamed_fasta)):
        return {}
    return {o: n for o, n in zip(_fasta_headers(orig_fasta), _fasta_headers(renamed_fasta))}


def _load_override_coords(cf, file_dict, gcols):
    """{(LOC, genome) -> (ref_chr_seqid, start, end)} = the EXACT own-genome annotation coordinates
    of the copy that genome CONTRIBUTED to a consensus locus, for placing Present:annot features
    precisely (instead of the unreliable flank-liftback gap). Built by joining the step-11 collapse
    dictionary (annots_incremental_dictionary.tsv: LOC -> original per-genome TE ids like
    'ID=<idx>_<genometag>_<superfam>#<family>|<len>') to that genome's own-coordinate annotation bed
    (slot 8, whose 4th-column leading token is the same <idx>), then remapping the own scaffold id
    into reference-chromosome space. The collapse recorded exactly which own copy contributed, so this
    also resolves the paralog ambiguity (flanks landing on different copies). A family-name check
    guards the index join. Empty {} when the dictionary is absent (e.g. groupby merge)."""
    dict_path = os.path.join(cf, "annots_incremental_dictionary.tsv")
    if not os.path.exists(dict_path):
        return {}
    tag2g = {g.replace("_", ""): g for g in gcols}       # species_tag (no "_") -> genome name
    own, seqmap = {}, {}                                  # per-genome: idx -> (seqid,s,e,idstr); own->ref map
    for g in gcols:
        slots = file_dict.get(g, [])
        bed = slots[8] if len(slots) > 8 else None
        d = {}
        if bed and os.path.exists(bed):
            with open(bed) as fh:
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) >= 4 and c[1].isdigit() and c[2].isdigit():
                        d[c[3].split("_", 1)[0]] = (c[0], int(c[1]), int(c[2]), c[3])
        own[g] = d
        seqmap[g] = _build_seqmap(slots[0] if len(slots) > 0 else None,
                                  slots[11] if len(slots) > 11 else None)
    out = {}
    with open(dict_path) as fh:
        next(fh, None)  # header
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            loc = c[3].split("|")[0]
            for oid in c[5].split(","):
                oid = oid.strip()
                if oid.startswith("ID="):
                    oid = oid[3:]
                parts = oid.split("_")
                if len(parts) < 2:
                    continue
                g = tag2g.get(parts[1])
                if g is None or (loc, g) in out:          # first contributing copy per (loc, genome)
                    continue
                rec = own[g].get(parts[0])
                fam = oid.split("#")[-1].split("|")[0] if "#" in oid else ""
                if rec is None or (fam and fam not in rec[3]):   # index join + family sanity check
                    continue
                sq, s, e, _ = rec
                out[(loc, g)] = (seqmap[g].get(sq, sq), s, e)
    return out


def _collect_burrite_features(matrix_path, transfer2_dir, file_dict=None, cf=None):
    """Per-genome BurriTE consensus features (Present/Absent/Ambiguous) from the matrix + each
    genome's flank-liftback coords. Returns (gcols, feats_by_g, skipped_by_g); each feature dict:
    {seqid,start,end,ftype,fid,family,superfamily,tags,coord_source} where tags = {Present,Absent,
    Ambiguous,Not_mapped} -> comma-joined OTHER-genome names ("none" if empty). ftype encodes the
    focal state: Present->transposable_element, Absent->insertion_site. Ambiguous (and Not_mapped)
    loci are NOT emitted as features (though other genomes' Ambiguous states still populate the tag).
    Coordinates (coord_source): confident genotypes use the real inter-flank span start=F1_EndR,
    end=F2_StartR ("flank_gap"). Present:annot overrides (the flank-liftback gap is unreliable there
    -- Ambiguous/NA/Absent, often a nonsensical multi-Mb span from flanks landing on different
    paralogs) instead take the EXACT own-genome annotation coordinates of the copy that genome
    contributed, looked up via the collapse dictionary + own-coord bed ("annotation", needs
    file_dict + cf); if that lookup fails they fall back to rebuilding from a lifted flank + the
    known TE length ("annot_anchor": F1_End..F1_End+L or F2_Start-L..F2_Start). Not_mapped loci are
    dropped; a locus that is neither dictionary-placeable nor flank-placeable is counted skipped."""
    def norm(s):
        s = str(s)
        if s.startswith("Present"):
            return "Present"
        if s.startswith("NA"):
            return "Not_mapped"
        return s if s in ("Absent", "Ambiguous") else "Not_mapped"

    # Ambiguous loci are excluded from the per-genome GFFs (only Present + Absent are emitted);
    # a state not in this map (Ambiguous / Not_mapped) is dropped.
    ftype = {"Present": "transposable_element", "Absent": "insertion_site"}

    def num(x):
        return int(x) if x.lstrip("-").isdigit() else None

    with open(matrix_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Incremental matrices have seven fixed columns, whereas the legacy
        # groupby matrix has six.  The Python v3 driver always sliced at seven,
        # silently dropping the first genome in groupby mode.  Nextflow v4
        # supports both layouts explicitly.
        fixed_cols = 7 if len(header) > 3 and header[3] == "LOC_ID" else 6
        gcols = header[fixed_cols:]
        rows = []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) <= fixed_cols:
                continue
            raw = {g: c[fixed_cols + i] for i, g in enumerate(gcols)}
            states = {g: norm(v) for g, v in raw.items()}
            if fixed_cols == 7:
                loc, fam, sup = c[3], c[4], c[5]
            else:
                loc, fam, sup = c[3], "NA", c[4]
            rows.append((loc, fam, sup, states, raw))

    # per-genome flank-liftback coords: LOC -> (chrom, F1_EndR, F2_StartR, TE_Length); a coord is
    # None when that flank did not lift (kept, so single-flank Present:annot loci can still be placed).
    coords = {}
    for g in gcols:
        d = {}
        path = os.path.join(transfer2_dir, f"{g}_Transfers_TE.txt")
        if os.path.exists(path):
            with open(path) as fh:
                next(fh, None)  # header
                for line in fh:
                    t = line.rstrip("\n").split("\t")
                    if len(t) < 7 or t[0] == "NA":               # both flanks failed -> unplaceable
                        continue
                    f1e, f2s = num(t[2]), num(t[3])
                    if f1e is None and f2s is None:
                        continue
                    coord = (t[0], f1e, f2s, num(t[6]))
                    d[t[5]] = coord
                    # Legacy groupby transfer tables identify loci as
                    # FM<n>_<superfamily>, while their matrix stores FM<n> and
                    # superfamily in separate columns.  Index both forms so
                    # per-genome GFF generation can recover the coordinates.
                    if t[5].startswith("FM") and "_" in t[5]:
                        d.setdefault(t[5].split("_", 1)[0], coord)
        coords[g] = d

    # exact own-annotation coords for Present:annot loci (dictionary join); {} if unavailable
    override_coords = _load_override_coords(cf, file_dict, gcols) if (file_dict and cf) else {}

    feats_by_g, skipped_by_g = {}, {}
    for g in gcols:
        others = [x for x in gcols if x != g]
        feats, skipped = [], 0
        for loc, fam, sup, states, raw in rows:
            st = states.get(g)
            if st not in ftype:                           # Not_mapped (or missing) -> skip
                continue
            override = raw[g].startswith("Present:")
            oc = override_coords.get((loc, g)) if override else None
            if oc is not None:
                # Present:annot: use the EXACT own-annotation coords of the copy this genome
                # contributed (from the collapse dictionary) -- resolves the paralog ambiguity and
                # is bp-exact, unlike the unreliable flank-liftback gap.
                seqid, s, e = oc
                coord_source = "annotation"
            else:
                cd = coords[g].get(loc)
                if cd is None:                            # no dict coords and both flanks NA
                    skipped += 1
                    continue
                seqid, f1e, f2s, L = cd
                if override:
                    # fallback: rebuild from a lifted flank + known TE length (flanks abut the TE)
                    Ln = L if L is not None else 0
                    s, e = (f1e, f1e + Ln) if f1e is not None else (f2s - Ln, f2s)
                    coord_source = "annot_anchor"
                elif f1e is not None and f2s is not None:
                    s, e = f1e, f2s                       # confident genotype: real inter-flank span
                    coord_source = "flank_gap"
                else:
                    skipped += 1                          # non-override needs both flanks (shouldn't happen)
                    continue
            if e < s:
                s, e = e, s
            s = max(1, s)
            buckets = {"Present": [], "Absent": [], "Ambiguous": [], "Not_mapped": []}
            for o in others:
                buckets[states[o]].append(o)
            tags = {k: (",".join(v) if v else "none") for k, v in buckets.items()}
            feats.append({"seqid": seqid, "start": s, "end": e, "ftype": ftype[st],
                          "fid": loc, "family": fam, "superfamily": sup, "tags": tags,
                          "coord_source": coord_source, "origin": "BurriTE", "prec": 3})
        feats_by_g[g] = feats
        skipped_by_g[g] = skipped
    return gcols, feats_by_g, skipped_by_g


def write_per_genome_gffs(matrix_path, transfer2_dir, out_dir, file_dict=None, cf=None):
    """Step 16 (--final_annotation burrite): one GFF per genome listing every consensus TE that is
    Present or Absent in that genome (Ambiguous and Not_mapped excluded). Feature type encodes the
    focal state (transposable_element / insertion_site); attributes carry ID=LOC_ID, family,
    superfamily and four tags listing the OTHER genomes in each state ("none" if empty)."""
    gcols, feats_by_g, skipped_by_g = _collect_burrite_features(matrix_path, transfer2_dir,
                                                                file_dict, cf)
    os.makedirs(out_dir, exist_ok=True)
    for g in gcols:
        outp = os.path.join(out_dir, f"{g}_TEs.gff")
        n = {"transposable_element": 0, "insertion_site": 0}
        with open(outp, "w") as out:
            out.write("##gff-version 3\n")
            for f in feats_by_g[g]:
                t = f["tags"]
                attrs = (f"ID={f['fid']};family={f['family']};superfamily={f['superfamily']};"
                         f"coord_source={f['coord_source']};"
                         f"Present={t['Present']};Absent={t['Absent']};"
                         f"Ambiguous={t['Ambiguous']};Not_mapped={t['Not_mapped']}")
                out.write(f"{f['seqid']}\tBurriTE\t{f['ftype']}\t{f['start']}\t{f['end']}"
                          f"\t.\t.\t.\t{attrs}\n")
                n[f["ftype"]] += 1
        msg = (f"INFO: step 16: {g}: {n['transposable_element']} Present / {n['insertion_site']} "
               f"Absent (Ambiguous excluded) -> {os.path.basename(outp)}")
        if skipped_by_g[g]:
            msg += f" ({skipped_by_g[g]} skipped: no mapped flank coords)"
        print(msg)


def _parse_source_features(final_annot_path, seqmap, min_len, max_len):
    """RepeatMasker + GraffiTE features from a genome's <sp>_final_annot.gff, length-filtered like
    step 7 and remapped into reference-chromosome coordinates. Origin is read from the GFF source
    column: 'RepeatMasker' vs 'Liftoff' (=GraffiTE). Returns (rm_feats, graffite_feats)."""
    rm, gf = [], []
    if not os.path.exists(final_annot_path):
        return rm, gf
    idx = 0
    with open(final_annot_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            try:
                s, e = int(c[3]), int(c[4])
            except ValueError:
                continue
            if e < s:
                s, e = e, s
            if (e - s) < min_len or (max_len is not None and (e - s) > max_len):
                continue                                    # same filter as filter_te_table (step 7)
            attrs = c[8]
            feat = {"seqid": seqmap.get(c[0], c[0]), "start": s, "end": e,
                    "family": _gff_attr(attrs, "Name") or "NA",
                    "superfamily": _gff_attr(attrs, "Class") or "NA",
                    "fid": _gff_attr(attrs, "ID") or f"f{idx}"}
            idx += 1
            if c[1] == "RepeatMasker":
                feat["origin"] = "RepeatMasker"; feat["prec"] = 1; rm.append(feat)
            elif c[1] == "Liftoff":
                feat["origin"] = "GraffiTE"; feat["prec"] = 2; gf.append(feat)
    return rm, gf


def _reciprocal_overlap(a, b, frac):
    """True if a and b overlap by >= frac of BOTH lengths (reciprocal overlap)."""
    ov = min(a["end"], b["end"]) - max(a["start"], b["start"])
    if ov <= 0:
        return False
    la, lb = a["end"] - a["start"], b["end"] - b["start"]
    return la > 0 and lb > 0 and ov >= frac * la and ov >= frac * lb


def _dedup_features(feats, frac=0.5):
    """Drop a feature only when it reciprocally overlaps (>= frac both ways) a kept feature from a
    STRICTLY higher-precedence source, so the same TE is not reported twice ACROSS sources while
    each source keeps its own annotation granularity (same-source overlaps are never collapsed).
    Precedence: BurriTE(3) > GraffiTE(2) > RepeatMasker(1). Per seqid: cluster by transitive
    any-overlap (cheap candidate generation), then within each cluster process by precedence desc,
    length desc; a lower-precedence feature that matches a higher-precedence kept feature is
    absorbed (its origin added to that host's merged_sources)."""
    by_seq = defaultdict(list)
    for f in feats:
        by_seq[f["seqid"]].append(f)
    kept = []
    for fs in by_seq.values():
        fs.sort(key=lambda f: (f["start"], f["end"]))
        clusters, cur, cur_max = [], [], None
        for f in fs:
            if cur and f["start"] < cur_max:
                cur.append(f); cur_max = max(cur_max, f["end"])
            else:
                if cur:
                    clusters.append(cur)
                cur, cur_max = [f], f["end"]
        if cur:
            clusters.append(cur)
        for cl in clusters:
            cl.sort(key=lambda f: (-f["prec"], -(f["end"] - f["start"])))
            local = []
            for f in cl:
                # absorb only into an already-kept feature from a HIGHER-precedence source
                host = next((k for k in local
                             if k["prec"] > f["prec"] and _reciprocal_overlap(f, k, frac)), None)
                if host is None:
                    f["merged_sources"] = {f["origin"]}
                    local.append(f)
                else:
                    host["merged_sources"].add(f["origin"])
            kept.extend(local)
    return kept


def write_per_genome_all_gffs(matrix_path, transfer2_dir, out_dir, file_dict,
                              min_length, max_length, cf=None, overlap_frac=0.5):
    """Step 16 (--final_annotation all): per-genome GFF combining the BurriTE consensus features
    with that genome's own RepeatMasker and GraffiTE annotations (length-filtered like step 7,
    remapped to reference-chromosome coords). Duplicates (same TE seen by >1 source) are collapsed
    by reciprocal overlap_frac overlap with precedence BurriTE > GraffiTE > RepeatMasker; the winner records
    merged_sources. RepeatMasker/GraffiTE features carry the Present/Absent/Ambiguous/Not_mapped
    tags = not_applicable (those cross-genome states are only defined for BurriTE consensus loci)."""
    gcols, feats_by_g, skipped_by_g = _collect_burrite_features(matrix_path, transfer2_dir,
                                                                file_dict, cf)
    os.makedirs(out_dir, exist_ok=True)
    NA = "not_applicable"
    for g in gcols:
        slots = file_dict.get(g, [])
        seqmap = _build_seqmap(slots[0] if len(slots) > 0 else None,
                               slots[11] if len(slots) > 11 else None)
        final_annot = slots[7] if len(slots) > 7 else None
        rm, gft = _parse_source_features(final_annot, seqmap, min_length, max_length) \
            if final_annot else ([], [])

        kept = _dedup_features(feats_by_g[g] + gft + rm, frac=overlap_frac)
        kept.sort(key=lambda f: (f["seqid"], f["start"], f["end"]))

        outp = os.path.join(out_dir, f"{g}_all_TEs.gff")
        counts = defaultdict(int)
        with open(outp, "w") as out:
            out.write("##gff-version 3\n")
            for f in kept:
                if f["origin"] == "BurriTE":
                    t = f["tags"]
                    P, A, Am, Nm = t["Present"], t["Absent"], t["Ambiguous"], t["Not_mapped"]
                    ftype = f["ftype"]
                    cs = f["coord_source"]
                else:
                    P = A = Am = Nm = NA
                    ftype = "transposable_element"
                    cs = "annotation"                     # RM/GraffiTE keep their native coords
                ms = ",".join(sorted(f["merged_sources"]))
                attrs = (f"ID={f['fid']};family={f['family']};superfamily={f['superfamily']};"
                         f"source={f['origin']};merged_sources={ms};coord_source={cs};"
                         f"Present={P};Absent={A};Ambiguous={Am};Not_mapped={Nm}")
                out.write(f"{f['seqid']}\t{f['origin']}\t{ftype}\t{f['start']}\t{f['end']}"
                          f"\t.\t.\t.\t{attrs}\n")
                counts[f["origin"]] += 1
        print(f"INFO: step 16 (all): {g}: {sum(counts.values())} features "
              f"[BurriTE {counts['BurriTE']} / GraffiTE {counts['GraffiTE']} / "
              f"RepeatMasker {counts['RepeatMasker']}] -> {os.path.basename(outp)}"
              + (f" ({skipped_by_g[g]} BurriTE skipped: no flank coords)" if skipped_by_g[g] else ""))


def build_contributor_map(merge_method, cf, collapsed_bed):
    """{FM<n>_<superfamily> -> set(contributing genome names)} from the step-11 collapse
    membership, for EITHER merge method. Consumed by build_presence_absence_matrix to mark a
    genome Present:annot wherever its own annotation contributed the locus (a TE is present
    by direct annotation, more reliable than the reference->genome flank round-trip).

      - incremental: read directly from annots_incremental_collapsed.tsv (samples column).
      - groupby:     annots_groupby.bed carries the contributing-genome list (col 6), but the
                     FM index is only assigned later (collapse_and_filter_grouped_bed), so map
                     the collapsed loci back to it by (chrom,start,end,superfamily).
    """
    contrib = {}
    if merge_method == "incremental":
        # tsv: Chrm start end LOC_ID n_ind class len_min len_max superfamily family
        #      n_families families samples Nested  -> key = LOC<n>, samples = col 12.
        tsv = os.path.join(cf, "annots_incremental_collapsed.tsv")
        with open(tsv) as fh:
            next(fh, None)  # header
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 13:
                    contrib[c[3].split("|")[0]] = set(c[12].split(","))
    else:  # groupby
        grouped = os.path.join(cf, "annots_groupby.bed")
        coord2sp = {}
        with open(grouped) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 6:
                    coord2sp[(c[0], c[1], c[2], c[3])] = set(c[5].split(","))
        with open(collapsed_bed) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 4:
                    continue
                name = c[3].split("|")[0]                    # FM<n>_<superfamily>
                superfam = name.split("_", 1)[1] if "_" in name else name
                contrib[name] = coord2sp.get((c[0], c[1], c[2], superfam), set())
    return contrib


def build_locus_meta(cf):
    """{LOC<n> -> (family, superfamily)} from the incremental collapse tsv, for the matrix's
    family/superfamily columns (incremental only)."""
    meta = {}
    tsv = os.path.join(cf, "annots_incremental_collapsed.tsv")
    with open(tsv) as fh:
        next(fh, None)  # header: ... [8]=superfamily [9]=family ...
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 13:
                meta[c[3].split("|")[0]] = (c[9], c[8])       # (family, superfamily)
    return meta


if __name__ == '__main__':

    Installation_path = os.path.dirname(os.path.realpath(__file__))
    tools_path = Installation_path + "/tools/"

    print("\n#########################################################################")
    print("#                                                                       #")
    print("#                             BurriTE:                                  #")
    print("#            A wrapper pipeline to annotate and analyze TEs in          #")
    print("#                        a pangenomic context.                          #")
    print("#                                                                       #")
    print("#########################################################################\n")

    ### read parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assemblies', required=True, dest='assemblies',
                        help='List of alternative assemblies (to be used by GraffiTE)')
    parser.add_argument('-t', '--te_lib', required=True, dest='te_lib',
                        help='Path to the Transposable Element library. Required*')
    parser.add_argument('-r', '--reference', required=True, dest='reference',
                        help='Reference genome to be used by GraffiTE. Required*')
    parser.add_argument('-o', '--output', required=False, dest='outputdir',
                        help='Path to the output directory.')
    parser.add_argument('-c', required=False, dest='cores', default=-1,
                        help='cores to execute some steps in parallel')
    parser.add_argument('--graffite_image', '--graffite-image', required=False,
                        dest='graffite_image', default=None,
                        help='Path to a previously downloaded GraffiTE Singularity/Apptainer '
                             'image (for example, graffite_latest.sif). When omitted, BurriTE '
                             'keeps GraffiTE\'s current default container behaviour.')
    parser.add_argument('--graffite_tmpdir', '--graffite-tmpdir', required=False,
                        dest='graffite_tmpdir', default=None,
                        help='Shared writable temporary directory for GraffiTE container tasks. '
                             'Default: <output>/01_GraffiTE/tmp.')
    parser.add_argument('-v', required=False, dest='verbose', default='N',
                        help='Verbose? [Y or N]. Default=N')
    parser.add_argument("--chromosome_list", required=True, help="File with chromosome names")
    parser.add_argument("--flank_size", required=False, type=int, default=500,
                        help="Number of bp to extract upstream and downstream of TE")
    parser.add_argument("--min_length", required=False, type=int, default=0,
                        help="Minimum length for keep a TE copy. Default=0 (no lower restriction).")
    parser.add_argument("--max_length", required=False, type=int, default=None,
                        help="Maximum length for keep a TE copy. Default=unset (no upper restriction).")
    parser.add_argument("--merge_method", required=False, choices=["groupby", "incremental"],
                        default=None,
                        help="Step 11 consensus collapse. Default=incremental "
                             "(superfamily+overlap+length clustering with shared/singleton "
                             "labels and per-locus membership). 'groupby' is the legacy "
                             "superfamily+adjacency collapse (coarser: over-merges distinct "
                             "copies then drops length-inconsistent groups; needs "
                             "bedtools>=2.31).")
    parser.add_argument("--include_reference", required=True, choices=["Y", "N"],
                        help="REQUIRED (no default). Integrate the reference genome's own TE "
                             "annotation into the consensus catalogue so reference(-private) loci "
                             "are genotyped across all queries, and add a reference genotype column "
                             "to the presence/absence matrix. [Y or N]. When Y and --merge_method is "
                             "unset, the incremental method is used.")
    parser.add_argument("--merge_by", required=False, choices=["family", "superfamily"],
                        default="superfamily",
                        help="[incremental] clustering granularity: 'superfamily' (default, "
                             "Class-level e.g. DNA_RC) or 'family' (Name-level e.g. ATREP13, "
                             "stricter, splits a superfamily locus into its families; better "
                             "for within-species). family/superfamily are written as separate "
                             "columns and loci get opaque LOC<n> ids. (groupby is superfamily-only.)")
    parser.add_argument("--final_annotation", required=False, choices=["burrite", "all"],
                        default="burrite",
                        help="Per-genome final GFF (step 16). 'burrite' (default): the consensus "
                             "TEs Present/Absent/Ambiguous in each genome (<g>_TEs.gff). 'all': "
                             "additionally merges that genome's own RepeatMasker + GraffiTE "
                             "annotations (length-filtered, deduped by reciprocal "
                             "--final_annotation_dedup_overlap overlap; precedence "
                             "BurriTE>GraffiTE>RepeatMasker) into <g>_all_TEs.gff.")
    parser.add_argument("--final_annotation_dedup_overlap", required=False, type=float, default=0.5,
                        help="Step 16 (--final_annotation all): reciprocal-overlap fraction (0-1) "
                             "above which two features from different sources are treated as the "
                             "same TE and collapsed by precedence. Default=0.5 (50%%).")
    parser.add_argument("--lenthr", required=False, type=int, default=100,
                        help="[incremental] max internal length spread within a cluster (bp). Default=100")
    parser.add_argument("--minlen", required=False, type=int, default=0,
                        help="[incremental] drop input intervals shorter than this (bp). Default=0")
    parser.add_argument("--bedmin", required=False, type=int, default=100,
                        help="[incremental] min consensus median length kept for the transfer BED (bp). Default=100")

    options = parser.parse_args()
    assemblies = options.assemblies
    te_lib = options.te_lib
    reference = options.reference
    outputdir = options.outputdir
    cores = options.cores
    graffite_image = options.graffite_image
    graffite_tmpdir = options.graffite_tmpdir
    verbose = options.verbose
    min_length = int(options.min_length)
    # max_length unset (None) means no upper length cap.
    max_length = int(options.max_length) if options.max_length is not None else None
    flank_size = int(options.flank_size)
    chromosome_list = options.chromosome_list
    merge_method = options.merge_method
    merge_by = options.merge_by
    final_annotation = options.final_annotation
    final_overlap = float(options.final_annotation_dedup_overlap)
    lenthr = int(options.lenthr)
    minlen = int(options.minlen)
    bedmin = int(options.bedmin)

    if graffite_image is not None:
        graffite_image = os.path.abspath(os.path.expanduser(graffite_image))
        if not os.path.isfile(graffite_image):
            parser.error("The GraffiTE image does not exist or is not a regular file: "
                         + graffite_image)

    # --include_reference (Y/N) toggles integrating the reference's own TE annotation.
    include_reference = str(options.include_reference).upper() in ["Y", "YES"]
    # Incremental is the global default: it clusters by superfamily + overlap + length so
    # distinct copies are kept as separate loci (groupby merges by superfamily + positional
    # adjacency and can over-merge then drop length-inconsistent groups). groupby is opt-in
    # legacy (pass --merge_method groupby explicitly).
    if merge_method is None:
        merge_method = "incremental"

    # Name the reference the SAME way the ALT samples are named: their labels come from
    # assemblies.csv column 2 (e.g. "100059"), which is also how parse_graffite_to_repeatmasker
    # derives them from the VCF (leading token before the first "."). Using splitext here left
    # the reference as "100002.Chr_scaffolds" while the ALTs were "100059", so it leaked a
    # ".Chr_scaffolds" suffix into the matrix column, reference bed and RM_<ref> dir. Take the
    # leading dot-token so "100002.Chr_scaffolds.fa" -> "100002" and the format matches for all
    # genomes. (The pipeline already assumes genome names contain no "."; see the Step 0 TODO.)
    reference_name = os.path.basename(reference).split(".")[0]

    if outputdir is None:
        outputdir = "BurriTE_output"

    if verbose is None:
        verbose = False
    elif verbose.upper() not in ["Y", "N", "YES", "NO"]:
        print("ERROR: Verbose should be Y, N, YES, or NO. Found "+verbose+" instead.")
    elif verbose.upper() in ["Y", "YES"]:
        verbose = True
    elif verbose.upper() in ["N", "NO"]:
        verbose = False

    if cores is None or cores == -1:
        cores = int(psutil.cpu_count())
        print("MESSAGE: Missing threads parameter, using by default: " + str(cores))
    else:
        cores = int(cores)
        
    # Creation of a Dict with assemblies names and files, where:
    # 0: assembly, 1: repeatMasker.out, 2: OneCode.out, 3: GraffiTE.out, 4: OneCode.out.gff, 5: GraffiTE.out.gff,
    # 6: transfered GraffiTEs.gff, 7: final joined gffs, 8: BED filtered files, 9: Anchors files, 10: mapping files
    # 11: renamed assemblies (with chr names from reference), 12: Anchors files with ref chr names
    file_dict = {}
    assemblies_lines = [x.replace("\n", "") for x in open(assemblies, 'r').readlines()]
    for i in range(1, len(assemblies_lines)):
        line = assemblies_lines[i]
        file_dict[line.split(",")[1].replace("\n", "")] = [line.split(",")[0]]
    # adding the reference (key must match reference_name derived above, otherwise the
    # DEL SVs that parse_graffite_to_repeatmasker assigns to reference_name have no dict slot)
    file_dict[reference_name] = [reference]

    ####################################################################################################################
    # Step 0: Checking
    ####################################################################################################################
    # TODO: To check that the assemblies dont have "." in their names

    ####################################################################################################################
    # Step 1: GraffiTE
    ####################################################################################################################
    create_output_folders(outputdir)
    create_output_folders(outputdir+"/01_GraffiTE")
    svf_file = run_GraffiTE(assemblies, te_lib, reference, cores,
                            outputdir+"/01_GraffiTE", verbose,
                            graffite_image, graffite_tmpdir)

    ####################################################################################################################
    # step 2. Run RepeatMasker for every species
    ####################################################################################################################
    create_output_folders(outputdir + "/02_RepeatMasker")
    file_dict = run_RepeatMasker(file_dict, te_lib, outputdir+"/02_RepeatMasker", cores, verbose)

    ####################################################################################################################
    # step 3. Run OneCodeToFindThemAll for every species
    ####################################################################################################################
    file_dict = run_OneCode(file_dict, outputdir+"/02_RepeatMasker", tools_path, verbose)

    ####################################################################################################################
    # step 4. Convert GraffiTE output and Repeatmasker to gff
    ####################################################################################################################
    file_dict = parse_graffite_to_repeatmasker(file_dict, reference_name, svf_file, outputdir+"/01_GraffiTE/3_TSD_search/")
    create_output_folders(outputdir + "/03_annotation_gff")
    file_dict = RM_out_2_gff(file_dict, outputdir+"/03_annotation_gff/")

    ####################################################################################################################
    # step 5. Transfer GraffiTE coordinates from reference to assemblies
    ####################################################################################################################
    file_dict = transfer_graffite_to_assemblies(file_dict, reference_name, outputdir + "/04_graffite_liftover/", verbose)

    ####################################################################################################################
    # step 6. Join GraffiTE and RepeatMasker annotations (without overlapping)
    file_dict = join_TE_annotations(file_dict, outputdir + "/05_joined_annotation/")
    ####################################################################################################################

    ####################################################################################################################
    # step 7. Prepare BED file
    create_output_folders(outputdir + "/06_lift_to_ref")
    file_dict = filter_te_table(file_dict, outputdir + "/06_lift_to_ref", min_length, max_length)
    ####################################################################################################################

    ####################################################################################################################
    # Step 8: Generate flanking regions and GFF
    file_dict = generate_flanking_regions(file_dict, flank_size, outputdir+ "/06_lift_to_ref")
    ####################################################################################################################

    ####################################################################################################################
    # Step 9: Separate genomes by corresponding chromosomes
    file_dict = run_minimap(file_dict, reference_name, outputdir+ "/06_lift_to_ref", cores)
    file_dict = parse_paf_unique_assignments(file_dict, reference_name, outputdir+ "/06_lift_to_ref")

    ####################################################################################################################
    # Step 10: Map and transfer anchors to the reference genome
    print(10)
    bed_file_folder = transfer_anchors(chromosome_list, file_dict, reference_name, cores, verbose, outputdir+ "/06_lift_to_ref")

    ####################################################################################################################
    # Step 11: Collapse the per-assembly calls into one consensus locus catalogue.
    # Two selectable methods; both write a 4-col FM<n>_<fam>|<len> BED (collapsed_output)
    # that the rest of the pipeline (steps 12-15) consumes unchanged.
    print(11)
    cf = outputdir + "/07_merge"
    create_output_folders(cf)
    if merge_method == "groupby" and merge_by == "family":
        print("WARNING: --merge_by family is only implemented for --merge_method incremental; "
              "groupby will collapse at superfamily level.")
    print(f"INFO: step 11 merge_method = {merge_method} | merge_by = {merge_by} | "
          f"include_reference = {include_reference}")

    # --include_reference: reformat the reference's own TE annotation into a per-genome BED
    # in bed_files/ so it is merged/collapsed together with the ALT-lifted loci below.
    ref_bed = None
    if include_reference:
        ref_bed = build_reference_bed(file_dict, reference_name, min_length, max_length,
                                      bed_file_folder)

    if merge_method == "incremental":
        sample_names = [a for a in file_dict.keys() if a != reference_name]
        bed_paths = [f"{bed_file_folder}/{a}.bed" for a in sample_names]
        if include_reference:
            # add the reference as its own individual (its private loci become singletons)
            sample_names = [reference_name] + sample_names
            bed_paths = [ref_bed] + bed_paths
        collapsed_output = incremental_collapse(bed_paths, sample_names, cf,
                                                "annots_incremental", lenthr, minlen, bedmin,
                                                merge_by=merge_by)
    else:  # "groupby" (validated); reference bed (if any) is picked up via bed_files/*.bed
        merged_bed = cf + "/annots_merged.bed"
        grouped_bed = cf + "/annots_groupby.bed"
        collapsed_output = cf + "/annots_groupby_collapsed.bed"
        merge_beds_to_manaus_nhv(bed_file_folder, merged_bed)
        run_bedtools_groupby(merged_bed, grouped_bed)
        collapse_and_filter_grouped_bed(grouped_bed, collapsed_output)

    # Collapse membership {locus -> contributing genomes} for the annotation-Present override,
    # and (incremental) the {LOC -> (family, superfamily)} metadata for the matrix columns.
    contrib = build_contributor_map(merge_method, cf, collapsed_output)
    loc_meta = build_locus_meta(cf) if merge_method == "incremental" else None

    ####################################################################################################################
    print(12)
    # Step 12: Extract flanking gffs and chrs from original assemblies
    for assembly in file_dict.keys():
        if assembly != reference_name:
            with open(chromosome_list) as f:
                chrom_list = [line.strip() for line in f]
            generate_flanking_gffs(collapsed_output, outputdir + "/07_merge", flank_size, chrom_list)
            # Use the chr-renamed assembly (slot 11), whose sequences carry the reference
            # chromosome names (e.g. CC1.8.Chr01); slot 0 is the original assembly whose
            # scaffolds have their own names, so nothing matched chrom_list and no FASTA
            # was written -> step 13 had no target sequence.
            extract_chromosome_fasta(file_dict[assembly][11], chrom_list, outputdir + "/07_merge", assembly)

    ####################################################################################################################
    print(13)
    # Step 13: Transfer the consensus FM flanks back from the reference onto each assembly
    # (the reverse of step 10) so every assembly is genotyped at the SAME unified loci.
    #   liftoff target    = the assembly chromosome (07_merge/<assembly>_<chrm>.fasta)
    #   liftoff reference  = the reference genome (annotation coords)
    #   annotation (-g)    = the consensus flanks  07_merge/FM_<chrm>.gff
    #   -f features        = 07_merge/<chrm>_gff_features.txt (type "repeat_element",
    #                        which matches FM_<chrm>.gff; the 04_graffite_liftover one is "repeat_region")
    # Output (per assembly x chrm) -> 07_merge/<assembly>_<chrm>_liffout.gff[_polished].
    step13_jobs = []
    per_job_threads = max(1, int(cores) // LIFTOFF_MAX_PARALLEL)
    for assembly in file_dict.keys():
        if assembly != reference_name:
            for chrm in chrom_list:
                fm_gff   = outputdir + f"/07_merge/FM_{chrm}.gff"
                query_fa = outputdir + f"/07_merge/{assembly}_{chrm}.fasta"
                features = outputdir + f"/07_merge/{chrm}_gff_features.txt"
                missing = next((p for p in (fm_gff, query_fa, features) if not os.path.exists(p)), None)
                if missing:
                    print(f"WARNING: Step 13 skipped for {assembly} {chrm}: missing {missing}.")
                    continue
                # Pre-index the (shared) reference genome and per-job query fasta serially so the
                # concurrent liftoff jobs never race building the .fai (see ensure_faidx).
                ensure_faidx(reference)
                ensure_faidx(query_fa)
                step13_jobs.append((fm_gff, query_fa, reference, features,
                                    f"{outputdir}/07_merge/{assembly}_{chrm}",
                                    outputdir + "/07_merge", per_job_threads, verbose))
    if step13_jobs:
        print(f"INFO: running {len(step13_jobs)} step-13 liftoff jobs, up to {LIFTOFF_MAX_PARALLEL} "
              f"in parallel ({per_job_threads} internal threads each)")
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(LIFTOFF_MAX_PARALLEL, len(step13_jobs))) as ex:
            list(ex.map(lambda a: run_liftoff(*a), step13_jobs))

    ####################################################################################################################
    print(14)
    # Step 14: Evaluate -- genotype each consensus FM locus in each assembly from the
    # step-13 lifted flanks. Concatenate this assembly's per-chromosome lifted polished
    # GFFs into one file, then evaluate it against the consensus FM list (collapsed_output).
    # 08_liftback holds the step 13/14 liftback intermediates (the concatenated
    # per-assembly _liffout_polished.gff and the _Transfers_TE.txt genotype tables); the
    # final deliverables (matrix + per-genome GFFs) stay in 09_final.
    transfer2_dir = outputdir + "/08_liftback"
    create_output_folders(transfer2_dir)
    create_output_folders(outputdir + "/09_final")
    for assembly in file_dict.keys():
        if assembly != reference_name:
            parts = [outputdir + f"/07_merge/{assembly}_{chrm}_liffout.gff_polished"
                     for chrm in chrom_list]
            parts = [p for p in parts if os.path.exists(p)]
            if not parts:
                print(f"WARNING: Step 14 skipped for {assembly}: no lifted polished GFF from step 13.")
                continue
            polished_gff_path = transfer2_dir + f"/{assembly}_liffout_polished.gff"
            with open(polished_gff_path, "w") as out:
                for p in parts:
                    with open(p) as fh:
                        shutil.copyfileobj(fh, out)
            final_output = transfer2_dir + f"/{assembly}_Transfers_TE.txt"
            evaluate_te_transfer(collapsed_output, polished_gff_path, final_output)

    ####################################################################################################################
    print(15)
    # Step 15: Join the per-assembly genotypes into one presence/absence matrix
    # (one row per consensus FM locus, one column per assembly).
    matrix_assemblies = [a for a in file_dict.keys() if a != reference_name]

    # --include_reference: genotype the reference itself (gap = FM locus end - start, since
    # its flanks are native) and add it as the first matrix column.
    if include_reference:
        ref_transfers = transfer2_dir + f"/{reference_name}_Transfers_TE.txt"
        evaluate_reference_transfer(collapsed_output, ref_transfers, flank_size)
        matrix_assemblies = [reference_name] + matrix_assemblies

    matrix_path = outputdir + "/09_final/Presence_Absence_matrix.txt"
    build_presence_absence_matrix(collapsed_output, matrix_assemblies,
                                  transfer2_dir, matrix_path,
                                  contrib=contrib, loc_meta=loc_meta)
    print(f"INFO: presence/absence matrix written to {matrix_path}")

    ####################################################################################################################
    print(16)
    # Step 16: per-genome GFF of the TEs Present/Absent in each genome (Ambiguous excluded; coords
    # from the collapse dictionary for annotation overrides else the flank liftback; feature type
    # encodes the focal state: transposable_element / insertion_site), tagging the other genomes as
    # Present/Absent/Ambiguous/Not_mapped. --final_annotation all also folds in that genome's own
    # RepeatMasker + GraffiTE annotations (deduped, precedence BurriTE>GraffiTE>RepeatMasker).
    per_genome_dir = outputdir + "/09_final/per_genome_gff"
    if final_annotation == "all":
        write_per_genome_all_gffs(matrix_path, transfer2_dir, per_genome_dir, file_dict,
                                  min_length, max_length, cf, overlap_frac=final_overlap)
    else:
        write_per_genome_gffs(matrix_path, transfer2_dir, per_genome_dir, file_dict, cf)
    print(f"INFO: per-genome TE GFFs ({final_annotation}) written to {per_genome_dir}")