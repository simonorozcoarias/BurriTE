# BurriTE
Scalable pan-annotation of transposable elements for population genomics

**BurriTE** builds transposable-element (TE) pan-annotations from a reference genome and multiple alternative assemblies. It integrates polymorphic TE calls from GraffiTE with homology-based annotations from RepeatMasker and One Code To Find Them All, reconciles orthologous TE loci through two rounds of flank transfer, and produces presence/absence matrices and enriched per-genome GFF3 files.

## Main features

- Runs RepeatMasker and One Code independently for every assembly and reference.
- Integrates polymorphic insertions detected by GraffiTE.
- Transfers TE-flanking anchors between assemblies and the reference with Minimap2 and Liftoff.
- Assigns persistent pan-locus identifiers across all genomes.
- Supports family- or superfamily-level locus clustering.
- Distinguishes `Present`, `Absent`, `Ambiguous` and `NA` genotypes.
- Generates a multi-genome presence/absence matrix.
- Produces enriched GFF3 annotations for every genome.
- Validates per-genome BED and transfer tables before global aggregation.
- Supports an existing GraffiTE VCF or a previously downloaded Apptainer/SIF image.
- Supports task-level caching and restart with Nextflow `-resume`.

## Workflow

For `N` alternative assemblies, one reference and `C` chromosomes, BurriTE can
schedule the following tasks:

| Stage | Number of tasks | Unit of parallelism |
|---|---:|---|
| GraffiTE | 1 nested workflow | complete pangenome |
| RepeatMasker | `N + 1` | genome |
| One Code To Find Them All | `N + 1` | genome |
| Annotation conversion and filtering | `N + 1` | genome |
| Minimap2 assembly alignment | `N` | alternative assembly |
| First Liftoff round | `N × C` | assembly × chromosome |
| Per-genome BED construction | `N` | alternative assembly |
| Pan-locus collapse | 1 | complete pangenome |
| Second Liftoff round | up to `N × C` | assembly × chromosome |
| Genotype evaluation | `N`, plus reference when enabled | genome |
| Presence/absence matrix | 1 | complete pangenome |
| Final GFF3 | `N`, or `N + 1` | genome |

The global pan-locus collapse and matrix construction are synchronization barriers because they require results from every genome.

GraffiTE remains a nested Nextflow workflow. Under the `slurm` profile, its launcher runs alongside the main Nextflow controller and GraffiTE submits its own cluster jobs using `--graffite_profile` (default: `cluster`).

## Requirements

- Linux.
- Java compatible with Nextflow 23.10.1.
- Nextflow 23.10.1. This version is pinned because GraffiTE uses the legacy
  Nextflow syntax parser.
- Python 3.10.
- RepeatMasker 4.2.4.
- Liftoff 1.6.3.
- Minimap2.
- SAMtools.
- BEDTools 2.31 or newer.
- Perl.
- Apptainer or Singularity on every node used by GraffiTE.
- The BurriTE-compatible One Code To Find Them All scripts.

## Installation

Clone or download this repository and create the Conda environment:

```bash
cd BurriTE
mamba env create -f envs/burrite.yml
conda activate burrite
```

If `mamba` is unavailable:

```bash
conda env create -f envs/burrite.yml
conda activate burrite
```

Confirm the principal programs:

```bash
nextflow -version
RepeatMasker -version
liftoff -h | head
minimap2 --version
samtools --version | head
bedtools --version
```

The Nextflow version should be:

```text
N E X T F L O W  ~  version 23.10.1
```

### Download the GraffiTE container once

Downloading the image before running BurriTE is recommended on clusters:

```bash
mkdir -p apptainer/images
apptainer pull apptainer/images/graffite_latest.sif \
    docker://cgroza/graffite:latest
```

Verify that the image can be opened:

```bash
apptainer exec apptainer/images/graffite_latest.sif minimap2 --version
```

Store the image on a shared filesystem visible from every SLURM node and pass
its absolute path with `--graffite_image`. If this parameter is omitted,
GraffiTE keeps its own default container configuration.

## Input files

BurriTE requires four inputs:

1. A CSV containing the alternative assemblies.
2. A curated TE library in FASTA format.
3. A reference genome in FASTA format.
4. A list of reference chromosome identifiers.

### Assemblies CSV

The CSV must contain the assembly path in the first column and the sample name
in the second column:

```csv
assembly,sample
/shared/project/genomes/C_arabica.fasta,C_arabica
/shared/project/genomes/C_brevipes.fasta,C_brevipes
/shared/project/genomes/C_congensis.fasta,C_congensis
```

Important rules:

- Do not include the reference genome in this CSV.
- Sample names must be unique.
- Sample names may contain letters, numbers, `_` and `-`.
- Do not use dots or whitespace in sample names.
- Assembly paths can be absolute or relative to the launch directory or CSV
  location. Absolute paths are recommended on clusters.

The reference sample name is inferred from the reference FASTA filename before
the first dot. For example:

```text
C_canephora_chr1.fasta  ->  C_canephora_chr1
```

This inferred name must not also occur in the assemblies CSV.

### TE library

Provide the curated TE consensus library used by both GraffiTE and
RepeatMasker:

```text
>TE_00000001#CLASSI/LTR/GYPSY
ACGT...
>TE_00000002#CLASSII/TIR/CACTA
ACGT...
```

The `#classification` suffix is strongly recommended because BurriTE uses TEfamily and superfamily information during pan-locus construction.

### Reference genome

The reference should be a FASTA assembly accessible from all compute nodes. Its chromosome identifiers must agree exactly with the chromosome list.

### Chromosome list

Provide one reference sequence identifier per line, without FASTA `>` symbols:

```text
CC1.8.Chr01
CC1.8.Chr02
CC1.8.Chr03
```

Every listed identifier must exist in the reference FASTA. BurriTE uses the Minimap2 PAF alignment to associate alternative-assembly contigs with these
reference identifiers.

## Quick start

### SLURM with a local GraffiTE image

```bash
nextflow run main.nf \
    -profile slurm \
    --assemblies test_data/assemblies.csv \
    --te_lib test_data/MCH_TElib_raw_Arabica_sgC_filtered_NR_oneCode.fasta \
    --reference test_data/C_canephora_chr1.fasta \
    --chromosome_list test_data/chr_list.txt \
    --graffite_tmpdir /shared/projects/my_project/tmp/graffite \
    --graffite_image /shared/projects/my_project/containers/graffite_latest.sif \
    --include_reference Y \
    --repeatmasker_cpus 40 \
    --minimap_cpus 40 \
    --liftoff_cpus 40 \
    --graffite_cpus 40 \
    --outdir BurriTE_output
```

Create the GraffiTE temporary directory beforehand:

```bash
mkdir -p /shared/projects/my_project/tmp/graffite
chmod 700 /shared/projects/my_project/tmp/graffite
```

The directory must be writable, have enough free space and be visible from all nodes used by GraffiTE. Avoid node-local `/tmp` unless the cluster guarantees that the same path exists inside every container task.

### Reuse a completed GraffiTE result

The most reliable way to skip GraffiTE completely is to provide its existing
VCF:

```bash
nextflow run main.nf \
    -profile slurm \
    --assemblies test_data/assemblies.csv \
    --te_lib test_data/MCH_TElib_raw_Arabica_sgC_filtered_NR_oneCode.fasta \
    --reference test_data/C_canephora_chr1.fasta \
    --chromosome_list test_data/chr_list.txt \
    --graffite_vcf "$(realpath BurriTE_output/01_GraffiTE/3_TSD_search/pangenome.vcf)" \
    --include_reference Y \
    --outdir BurriTE_output \
    -resume
```

This bypasses `RUN_GRAFFITE`; the supplied VCF becomes the input to `PARSE_GRAFFITE_VCF`.

### Local execution

For a small test or workstation, use an existing GraffiTE VCF:

```bash
nextflow run main.nf \
    -profile local \
    --assemblies examples/assemblies.csv \
    --te_lib /path/to/TE_library.fasta \
    --reference /path/to/reference.fasta \
    --chromosome_list examples/chromosomes.txt \
    --graffite_vcf /path/to/pangenome.vcf \
    --outdir BurriTE_test \
    -resume
```

The default GraffiTE profile targets a cluster, so `--graffite_vcf` is recommended for local tests.

## Running the Nextflow controller on a cluster

The main Nextflow process is a coordinator. It does not perform RepeatMasker, Minimap2 or Liftoff computations itself, so it does not need the sum of all task CPUs.

A typical controller allocation is:

- 1 CPU.
- 4–8 GB RAM.
- A long wall time covering the complete workflow.
- Access to the shared project, `work/`, output and GraffiTE temporary paths.

Run the controller in `screen`, `tmux` or a small long-running SLURM job, according to local cluster policy. The individual processes request their own resources and are submitted independently.

If the cluster automatically assigns memory per CPU and rejects explicit memory requests, remove or override the `memory` directives in `nextflow.config` for that installation. Do not increase the controller to 40 CPUs merely because the computational processes use 40 CPUs.

## Parameters

### Required parameters

| Parameter | Description |
|---|---|
| `--assemblies` | CSV containing `assembly,sample` for alternative assemblies |
| `--te_lib` | Curated TE-library FASTA |
| `--reference` | Reference-genome FASTA |
| `--chromosome_list` | Reference sequence IDs, one per line |

### GraffiTE parameters

| Parameter | Default | Description |
|---|---:|---|
| `--graffite_vcf` | unset | Existing `pangenome.vcf`; skips GraffiTE |
| `--graffite_image` | unset | Previously downloaded SIF/IMG container |
| `--graffite_tmpdir` | `<outdir>/01_GraffiTE/tmp` | Shared writable temporary directory |
| `--graffite_profile` | `cluster` | Profile passed to the nested GraffiTE workflow |
| `--graffite_cpus` | `40` | Cores passed to GraffiTE |

BurriTE currently runs GraffiTE with `graphaligner` and with read-based genotyping disabled (`--genotype false`).

### Analysis parameters

| Parameter | Default | Description |
|---|---:|---|
| `--include_reference` | `Y` | Include the reference annotation and reference-private loci |
| `--flank_size` | `500` | Flank length around each TE |
| `--min_length` | `0` | Minimum TE-copy length |
| `--max_length` | unset | Maximum TE-copy length |
| `--merge_method` | `incremental` | `incremental` or legacy `groupby` collapse |
| `--merge_by` | `superfamily` | Incremental clustering by `family` or `superfamily` |
| `--lenthr` | `100` | Maximum internal length spread during incremental clustering |
| `--minlen` | `0` | Minimum interval length accepted during collapse |
| `--bedmin` | `100` | Minimum median pan-locus length transferred back |
| `--final_annotation` | `burrite` | `burrite` or integrated `all` GFF mode |
| `--final_annotation_dedup_overlap` | `0.5` | Reciprocal-overlap threshold used by `all` mode |
| `--allow_empty_sample_beds` | `false` | Allow empty first-round BEDs; diagnostic/edge-case option |

`--allow_empty_sample_beds true` should only be used when an assembly genuinely contains no usable TE loci. It should not be used to bypass a Liftoff or GFF formatting error.

### Resource and SLURM parameters

| Parameter | Default | Description |
|---|---:|---|
| `--repeatmasker_cpus` | `16` | CPUs per RepeatMasker task |
| `--minimap_cpus` | `16` | CPUs per Minimap2 task |
| `--liftoff_cpus` | `4` | CPUs per Liftoff task |
| `--queue` | unset | SLURM partition/queue |
| `--cluster_options` | unset | Additional options passed to `sbatch` |

RMBlast uses approximately four threads per RepeatMasker `-pa` worker. BurriTE therefore derives `-pa` as `repeatmasker_cpus / 4` so the SLURM CPU request represents the approximate maximum thread usage.

Examples:

```bash
--queue normal
```

```bash
--cluster_options='--account=my_account --qos=normal'
```

## Output structure

```text
BurriTE_output/
├── 01_GraffiTE/
│   └── 3_TSD_search/pangenome.vcf
├── 02_RepeatMasker/
│   └── RM_<sample>/
├── 03_annotation_gff/
├── 04_graffite_liftover/
├── 05_joined_annotation/
├── 06_lift_to_ref/
│   ├── anchors/
│   ├── bed_files/
│   ├── filtered_annotation/
│   ├── liftoff/
│   ├── paf/
│   ├── renamed/
│   ├── transfer_<sample>/
│   └── validation/sample_bed_validation.tsv
├── 07_merge/
│   ├── collapsed.bed
│   └── liftback_liftoff/
├── 08_liftback/
│   ├── <sample>_Transfers_TE.txt
│   └── validation/transfer_validation.tsv
└── 09_final/
    ├── Presence_Absence_matrix.txt
    └── per_genome_gff/
        └── <sample>_TEs.gff
```

Nextflow also writes the following execution reports in the launch directory:

```text
BurriTE_v4_report.html
BurriTE_v4_timeline.html
BurriTE_v4_trace.txt
BurriTE_v4_dag.html
```

### Presence/absence matrix

`09_final/Presence_Absence_matrix.txt` contains one row per pan-locus and the following metadata columns:

```text
Chrm  start  end  LOC_ID  family  superfamily  TE_length  <genomes...>
```

Genome cells can contain:

| State | Meaning |
|---|---|
| `Present` | Inter-flank span is compatible with the TE length |
| `Present:annot` | Direct annotation confirms a contributing copy |
| `Absent` | Both flanks map with an empty or near-empty insertion interval |
| `Ambiguous` | Both flanks map, but the interval is not decisively present or absent |
| `NA:<reason>` | One or both flanks could not be evaluated reliably |

The default presence rule accepts an inter-flank span within 20% of the pan-locus TE length. An absence is called when the gap between flanks is 0–10 bp. Direct annotation evidence takes precedence over an uncertain flank round-trip and is recorded as `Present:annot`.

### Final GFF3 modes

With the default:

```bash
--final_annotation burrite
```

BurriTE writes `<sample>_TEs.gff`, containing pan-loci classified as present or absent for that genome. Present copies use the `transposable_element` feature type and absent sites use `insertion_site`.

With:

```bash
--final_annotation all
```

BurriTE writes `<sample>_all_TEs.gff` and additionally integrates the genome's own RepeatMasker and GraffiTE features. Overlapping calls are deduplicated with the precedence:

```text
BurriTE > GraffiTE > RepeatMasker
```

## Resume behavior

Restart a failed or interrupted execution with the same command and add:

```bash
-resume
```

For caching to work:

- Keep the same `work/` directory.
- Do not modify input files in place.
- Keep parameter values unchanged unless a rerun is intended.
- Avoid mixing `main.nf` and Python helper files from different BurriTE
  revisions.

Each genome and chromosome is cached independently. A failure in one Liftoff job does not require completed RepeatMasker or Minimap2 jobs to run again.

The nested GraffiTE launcher also uses `-resume`, with its work directory under `<outdir>/.graffite_work`. Nevertheless, supplying `--graffite_vcf` is the deterministic way to guarantee that an already completed GraffiTE analysis is
not relaunched.

## Reproducibility recommendations

- Keep Nextflow pinned to 23.10.1 for this BurriTE/GraffiTE combination.
- Record the BurriTE commit or release tag.
- Record the SHA-256 checksum of the GraffiTE image.
- Keep the input TE library and assemblies immutable during a run.
- Preserve `BurriTE_v4_trace.txt`, the execution report and validation tables.
- Use absolute shared paths for cluster inputs, containers and temporary data.

