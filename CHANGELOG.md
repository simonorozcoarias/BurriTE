# Changelog

## 4.0.0

- Reimplemented the BurriTE driver as a Nextflow DSL2 workflow.
- Parallelized RepeatMasker, OneCode and annotation preparation per genome.
- Parallelized both Liftoff rounds per alternative assembly and chromosome.
- Added independent CPU/memory/time labels and local/SLURM profiles.
- Preserved the default GraffiTE container behavior and the optional local
  `--graffite_image` and shared `--graffite_tmpdir` parameters.
- Added `--graffite_vcf` to reuse an existing GraffiTE result.
- Retained reference integration, incremental/groupby collapse, family or
  superfamily clustering, presence/absence matrices and both final-GFF modes.
- Fixed legacy groupby final-GFF parsing: its six fixed matrix columns and
  `FM<n>_<superfamily>` transfer identifiers are now handled correctly.
- Replaced the v3 OneCode logging path that referenced an undefined file handle
  with explicit per-task logs.
- Added fail-fast input validation, a JSON parameter schema, a stage smoke test
  and a complete Nextflow `-preview` topology test.
