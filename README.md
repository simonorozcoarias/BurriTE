# BurriTE
Pan-annotation of transposable elements for population 

## Parallel execution model

For `N` alternative assemblies, one reference and `C` chromosomes, the workflow
can create:

| Stage | Number of tasks | Unit of parallelism |
|---|---:|---|
| RepeatMasker | `N + 1` | assembly/reference |
| One Code To Find Them All | `N + 1` | assembly/reference |
| Annotation conversion/filtering | `N + 1` | assembly/reference |
| Minimap2 | `N` | alternative assembly |
| First Liftoff round | `N Ã— C` | assembly Ã— chromosome |
| Per-assembly BED construction | `N` | alternative assembly |
| Second Liftoff round | up to `N Ã— C` | assembly Ã— chromosome with consensus loci |
| Genotype evaluation | `N` | alternative assembly |
| Final GFF | `N`, or `N + 1` | genome |

The consensus collapse and presence/absence matrix are global barriers because
they require all per-genome results.

Version 4.0.1 fixes the assembly-by-chromosome channel expansion in both
Liftoff rounds. It also checks that no genome disappears silently before the
global collapse or final matrix is generated.

GraffiTE remains a nested Nextflow workflow. Its launcher runs locally on the
submission node and GraffiTE submits its own jobs using its `cluster` profile.
This avoids launching a second Nextflow controller inside a SLURM compute job.

## Requirements

- Linux cluster or workstation.
- Nextflow 23.10.1 and Java.
- RepeatMasker, Liftoff, Minimap2, samtools and bedtools >=2.31.
- Apptainer/Singularity available on every node used by GraffiTE.
- The BurriTE-compatible One Code To Find Them All scripts.

Create the portable environment with:

```bash
mamba env create -f envs/burrite.yml
conda activate burrite
```
## Input files

The assemblies CSV has two columns and does **not** include the reference:

```csv
assembly,sample
/shared/project/genome_A.fasta,genome_A
/shared/project/genome_B.fasta,genome_B
/shared/project/genome_C.fasta,genome_C
```

Sample IDs may contain letters, numbers, `_` and `-`, but not dots. FASTA paths
may be absolute or relative to the launch directory/CSV location.

The chromosome list contains one reference chromosome ID per line. Each ID must
exist in the reference FASTA and, after the Minimap2 renaming stage, in every
alternative assembly:

```text
Chr01
Chr02
Chr03
```

## SLURM execution

```bash
nextflow run main.nf \
    -profile slurm \
    --assemblies test_data/assemblies.csv \
    --te_lib test_data/MCH_TElib_raw_Arabica_sgC_filtered_NR_oneCode.fasta \
    --reference test_data/C_arabica_sgC_chr1.fasta \
    --chromosome_list test_data/chromosomes.txt \
    --graffite_image /shared/containers/graffite_latest.sif \
    --graffite_tmpdir /shared/scratch/$USER/graffite_tmp \
    --include_reference Y \
    --repeatmasker_cpus 40 \
    --minimap_cpus 40 \
    --liftoff_cpus 40 \
    --graffite_cpus 40 \
    --outdir BurriTE_output \
```

When `--graffite_image` is omitted, GraffiTE retains its original/default
container behavior. To reuse a completed GraffiTE run and skip the nested
workflow:

```bash
--graffite_vcf /path/to/01_GraffiTE/3_TSD_search/pangenome.vcf
```

## Local execution

For a small test or workstation:

```bash
nextflow run main.nf -profile local \
    --assemblies assemblies.csv \
    --te_lib library.fasta \
    --reference reference.fasta \
    --chromosome_list chromosomes.txt \
    --graffite_vcf existing_pangenome.vcf \
    --outdir BurriTE_test \
    -resume
```

Using `--graffite_vcf` is recommended for local tests because the default
GraffiTE profile is configured for a cluster.

## Important parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--include_reference` | `Y` | Include and genotype reference-private loci |
| `--allow_empty_sample_beds` | `false` | Permit an empty BED after first Liftoff; diagnostic use only |
| `--flank_size` | `500` | Flank length around each TE |
| `--min_length` | `0` | Minimum TE-copy length |
| `--max_length` | unset | Maximum TE-copy length |
| `--merge_method` | `incremental` | `incremental` or legacy `groupby` |
| `--merge_by` | `superfamily` | Incremental clustering by family/superfamily |
| `--lenthr` | `100` | Maximum internal length spread in a cluster |
| `--minlen` | `0` | Minimum input interval length during collapse |
| `--bedmin` | `100` | Minimum consensus median length transferred back |
| `--final_annotation` | `burrite` | `burrite` or integrated `all` GFF |
| `--final_annotation_dedup_overlap` | `0.5` | Reciprocal overlap for final source deduplication |
| `--repeatmasker_cpus` | `16` | SLURM CPUs per RepeatMasker task |
| `--minimap_cpus` | `16` | CPUs per Minimap2 task |
| `--liftoff_cpus` | `4` | CPUs per Liftoff task |
| `--graffite_cpus` | `40` | Cores passed to the nested GraffiTE workflow |
| `--queue` | unset | SLURM partition |
| `--cluster_options` | unset | Additional SLURM options |

RepeatMasker with RMBlast uses approximately four threads per `-pa` worker.
BurriTE therefore derives `-pa` as `repeatmasker_cpus / 4`; the requested SLURM
CPUs reflect the real maximum thread use.

## Output structure

```text
BurriTE_output/
      01_GraffiTE/
      02_RepeatMasker/
      03_annotation_gff/
      04_graffite_liftover/
      05_joined_annotation/
      06_lift_to_ref/
          bed_files/
          liftoff/
          validation/sample_bed_validation.tsv
      07_merge/
      08_liftback/
          validation/transfer_validation.tsv
      09_final/
          Presence_Absence_matrix.txt
          per_genome_gff/
```

Nextflow also writes `BurriTE_v4_report.html`, `BurriTE_v4_timeline.html`,
`BurriTE_v4_trace.txt` and `BurriTE_v4_dag.html` in the launch directory.

## Resume and failure behavior

Use `-resume` on every restart. Each assembly and chromosome is cached
independently, so a failed task does not force completed RepeatMasker, Minimap2
or Liftoff tasks to run again.

RepeatMasker, OneCode, Minimap2 and the first GraffiTE-coordinate Liftoff stage
fail fast. The two flank-transfer Liftoff stages preserve v3 behavior: a failed
mapping is represented by an empty polished GFF and becomes `NA`/unmappable in
the final matrix rather than terminating the entire pan-genome analysis.

Before collapse, BurriTE requires exactly one non-empty BED for every expected
genome and writes `06_lift_to_ref/validation/sample_bed_validation.tsv`. If a
BED is empty, inspect the corresponding first-round log in
`06_lift_to_ref/liftoff/`. Only use `--allow_empty_sample_beds true` when a
genome genuinely has no usable TE loci and that behavior is intended.

Before building the matrix, BurriTE checks that every sample has one transfer
table with exactly one row per collapsed locus. Genotype totals are written to
`08_liftback/validation/transfer_validation.tsv`; a missing or truncated table
stops the workflow instead of producing a silently empty final GFF.
ranges. A complete biological run still requires the programs listed under
Requirements and real assemblies/TE annotations.
