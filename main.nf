nextflow.enable.dsl = 2

/*
 * BurriTE v4
 * Pan-annotation of transposable elements with per-assembly and
 * per-assembly/chromosome parallelism.
 */

params.assemblies            = null
params.te_lib                = null
params.reference             = null
params.chromosome_list       = null
params.outdir                = 'BurriTE_output'
params.graffite_vcf          = null
params.graffite_image        = null
params.graffite_tmpdir       = null
params.graffite_profile      = 'cluster'
params.graffite_cpus         = 40
params.onecode_dir           = "${projectDir}/tools/OneCodeToFindThemAll"
params.include_reference     = 'Y'
params.flank_size            = 500
params.min_length            = 0
params.max_length            = null
params.merge_method          = 'incremental'
params.merge_by              = 'superfamily'
params.final_annotation      = 'burrite'
params.final_annotation_dedup_overlap = 0.5
params.allow_empty_sample_beds = false
params.lenthr                = 100
params.minlen                = 0
params.bedmin                = 100


def enabled(value) {
    value != null && ['1', 'TRUE', 'T', 'YES', 'Y'].contains(value.toString().toUpperCase())
}


def resolveCsvPath(String rawPath, java.nio.file.Path sheetParent) {
    def candidate = java.nio.file.Paths.get(rawPath).normalize()
    if( !candidate.isAbsolute() ) {
        def fromLaunch = java.nio.file.Paths.get(launchDir.toString()).resolve(candidate).normalize()
        def fromSheet = sheetParent.resolve(candidate).normalize()
        candidate = java.nio.file.Files.exists(fromLaunch) ? fromLaunch : fromSheet
    }
    file(candidate.toString(), checkIfExists: true)
}


process RUN_GRAFFITE {
    label 'orchestrator'
    tag 'GraffiTE'

    input:
    path assemblies_csv
    path te_lib
    path reference
    val output_root
    val temporary_dir
    val local_image

    output:
    path 'pangenome.vcf', emit: vcf
    path 'GraffiTE.log', emit: log

    script:
    def imageOption = local_image ? "-with-singularity '${local_image}'" : ''
    """
    set -euo pipefail

    python '${projectDir}/bin/burrite_stage.py' normalize-samplesheet \
        --input '${assemblies_csv}' \
        --base-dir '${launchDir}' \
        --output assemblies.absolute.csv

    mkdir -p '${output_root}/01_GraffiTE' '${output_root}/.graffite_work' '${temporary_dir}'
    chmod 700 '${temporary_dir}'

    bind_path="\${APPTAINER_BINDPATH:-}"
    if [[ -n "\${bind_path}" ]]; then
        bind_path="\${bind_path},${temporary_dir}"
    else
        bind_path='${temporary_dir}'
    fi

    env -u SINGULARITYENV_TMP -u SINGULARITYENV_TMPDIR \
        -u NXF_TASK_WORKDIR -u NXF_CHDIR \
        TMP='${temporary_dir}' \
        TMPDIR='${temporary_dir}' \
        APPTAINERENV_TMP='${temporary_dir}' \
        APPTAINERENV_TMPDIR='${temporary_dir}' \
        APPTAINER_BINDPATH="\${bind_path}" \
        NXF_SYNTAX_PARSER=v1 \
        nextflow run cgroza/GraffiTE \
            -profile '${params.graffite_profile}' \
            -work-dir '${output_root}/.graffite_work' \
            -resume \
            ${imageOption} \
            --assemblies assemblies.absolute.csv \
            --TE_library '${te_lib}' \
            --reference '${reference}' \
            --graph_method graphaligner \
            --genotype false \
            --cores ${params.graffite_cpus} \
            --out '${output_root}/01_GraffiTE' \
            > GraffiTE.log 2>&1

    test -s '${output_root}/01_GraffiTE/3_TSD_search/pangenome.vcf'
    cp '${output_root}/01_GraffiTE/3_TSD_search/pangenome.vcf' pangenome.vcf
    """
}


process PARSE_GRAFFITE_VCF {
    label 'small'
    tag 'GraffiTE VCF -> annotations'
    publishDir "${params.outdir}/01_GraffiTE/3_TSD_search", mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.tokenize('/')[-1] }

    input:
    path vcf
    val samples
    val reference_name

    output:
    path 'graffite_rm/*_GraffiTE.out', emit: records

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' graffite-vcf-to-rm \
        --vcf '${vcf}' \
        --samples '${samples}' \
        --reference-name '${reference_name}' \
        --output-dir graffite_rm
    """
}


process REPEATMASKER {
    label 'repeatmasker'
    tag "${sample}"
    publishDir "${params.outdir}/02_RepeatMasker", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fasta)
    path te_lib

    output:
    tuple val(sample), path("RM_${sample}"), emit: rm_dir

    script:
    // RMBlast uses four threads per RepeatMasker -pa worker.  Deriving -pa
    // from task.cpus keeps the SLURM allocation honest.
    def rmParallel = Math.max(1, (task.cpus as int).intdiv(4))
    """
    set -euo pipefail
    mkdir -p 'RM_${sample}'
    RepeatMasker \
        -famdb_dir '' \
        -pa ${rmParallel} \
        -dir 'RM_${sample}' \
        -lib '${te_lib}' \
        '${fasta}' \
        > 'RM_${sample}/RepeatMasker.stdout.log' \
        2> 'RM_${sample}/RepeatMasker.stderr.log'
    test -n "\$(find 'RM_${sample}' -maxdepth 1 -type f -name '*.out' -print -quit)"
    """
}


process ONECODE {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/02_RepeatMasker", mode: 'copy', overwrite: true,
        saveAs: { filename -> "RM_${sample}/${filename.tokenize('/')[-1]}" }

    input:
    tuple val(sample), path(rm_dir)
    path onecode_dir

    output:
    tuple val(sample),
          path('onecode/*.out.elem_sorted.csv.out'),
          path('onecode/dico.txt'),
          path('onecode/OneCode.log'), emit: result

    script:
    """
    set -euo pipefail
    test -f '${onecode_dir}/build_dictionary.pl'
    test -f '${onecode_dir}/one_code_to_find_them_all_but_sanely.pl'
    mkdir -p onecode
    rm_source=\$(find -L '${rm_dir}' -maxdepth 1 -type f -name '*.out' ! -name '*.elem*' | head -n 1)
    test -n "\${rm_source}"
    rm_out="onecode/\$(basename "\${rm_source}")"
    cp "\${rm_source}" "\${rm_out}"

    perl '${onecode_dir}/build_dictionary.pl' --rm "\${rm_out}" \
        > onecode/dico.txt 2> onecode/build_dictionary.stderr.log
    perl '${onecode_dir}/one_code_to_find_them_all_but_sanely.pl' \
        --rm "\${rm_out}" --ltr onecode/dico.txt --unknown --strict \
        > onecode/OneCode.log 2>&1
    test -f "\${rm_out}.elem_sorted.csv"
    awk '/^###/{sub(/^###/,""); print}' "\${rm_out}.elem_sorted.csv" \
        > "\${rm_out}.elem_sorted.csv.out"
    """
}


process CONVERT_ANNOTATIONS {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/03_annotation_gff", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(repeatmasker), path(graffite)

    output:
    tuple val(sample),
          path("${sample}_RepeatMasker.gff"),
          path("${sample}_GraffiTE.gff"), emit: gffs

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' annotation-gffs \
        --sample '${sample}' \
        --repeatmasker '${repeatmasker}' \
        --graffite '${graffite}'
    """
}


process LIFT_GRAFFITE {
    label 'liftoff'
    tag "${sample}"
    publishDir "${params.outdir}/04_graffite_liftover", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(assembly), path(graffite_gff)
    path reference

    output:
    tuple val(sample), path("${sample}_GraffiTE_lifted.gff"), emit: lifted

    script:
    """
    set -euo pipefail
    printf 'repeat_region\n' > features_file.txt
    mkdir -p liftoff_tmp
    if [[ -s '${graffite_gff}' ]]; then
        liftoff \
            -g '${graffite_gff}' \
            -o '${sample}_GraffiTE_lifted.gff' \
            -f features_file.txt \
            -dir liftoff_tmp \
            -p ${task.cpus} \
            '${assembly}' '${reference}' \
            > '${sample}_GraffiTE_liftoff.log' 2>&1
    else
        : > '${sample}_GraffiTE_lifted.gff'
    fi
    """
}


process COPY_REFERENCE_GRAFFITE {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/04_graffite_liftover", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(graffite_gff)

    output:
    tuple val(sample), path("${sample}_GraffiTE_lifted.gff"), emit: lifted

    script:
    """
    cp '${graffite_gff}' '${sample}_GraffiTE_lifted.gff'
    """
}


process JOIN_ANNOTATIONS {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/05_joined_annotation", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(repeatmasker_gff), path(graffite_gff)

    output:
    tuple val(sample), path("${sample}_final_annot.gff"), emit: annotation

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' join-annotations \
        --sample '${sample}' \
        --repeatmasker-gff '${repeatmasker_gff}' \
        --graffite-gff '${graffite_gff}'
    """
}


process FILTER_ANNOTATION {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/06_lift_to_ref/filtered_annotation", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(annotation)

    output:
    tuple val(sample), path("${sample}_final_annot.bed"), emit: bed

    script:
    def maxOption = params.max_length != null ? "--max-length ${params.max_length}" : ''
    """
    python '${projectDir}/bin/burrite_stage.py' filter-bed \
        --sample '${sample}' \
        --annotation '${annotation}' \
        --min-length ${params.min_length} \
        ${maxOption}
    """
}


process MAKE_ANCHORS {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/06_lift_to_ref/anchors", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(bed)

    output:
    tuple val(sample),
          path("${sample}_anchors${params.flank_size}.gff"),
          path("${sample}_anchor_id_map.tsv"), emit: anchors

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' make-anchors \
        --sample '${sample}' \
        --bed '${bed}' \
        --flank-size ${params.flank_size} \
        --gff-id-encoding short-map-v2
    """
}


process MINIMAP_ASSEMBLY {
    label 'minimap'
    tag "${sample}"
    publishDir "${params.outdir}/06_lift_to_ref/paf", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(assembly)
    path reference

    output:
    tuple val(sample), path("${sample}.paf"), emit: paf

    script:
    """
    minimap2 -x asm5 -t ${task.cpus} '${reference}' '${assembly}' \
        > '${sample}.paf' 2> '${sample}.minimap2.log'
    """
}


process RENAME_ASSEMBLY {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/06_lift_to_ref/renamed", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(assembly), path(anchors), path(anchor_map), path(paf)

    output:
    tuple val(sample),
          path("${sample}_chrRenamed.fasta"),
          path("${sample}_anchors_chrRenamed.gff"),
          path("${sample}_renamed_anchor_id_map.tsv"), emit: renamed

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' rename-by-paf \
        --sample '${sample}' \
        --fasta '${assembly}' \
        --anchors '${anchors}' \
        --paf '${paf}'
    cp '${anchor_map}' '${sample}_renamed_anchor_id_map.tsv'
    """
}


process EXTRACT_REFERENCE_CHROM {
    label 'small'
    tag "${chromosome}"

    input:
    val chromosome
    path reference

    output:
    tuple val(chromosome), path("reference_${chromosome}.fasta"), emit: chromosome_fasta

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' extract-reference-chrom \
        --reference '${reference}' \
        --chromosome '${chromosome}'
    """
}


process PREPARE_ALT_CHROM {
    label 'small'
    tag "${sample}:${chromosome}"

    input:
    tuple val(sample), val(chromosome), path(renamed_fasta), path(renamed_anchors),
          path(anchor_map)

    output:
    tuple val(sample), val(chromosome),
          path("${sample}_${chromosome}.fasta"),
          path("${sample}_${chromosome}_anchors.gff"),
          path("${sample}_${chromosome}_anchor_id_map.tsv"),
          path("${chromosome}_repeat_region.features"), emit: prepared

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' prepare-alt-chrom \
        --sample '${sample}' \
        --chromosome '${chromosome}' \
        --fasta '${renamed_fasta}' \
        --anchors '${renamed_anchors}'
    cp '${anchor_map}' '${sample}_${chromosome}_anchor_id_map.tsv'
    """
}


process LIFT_ANCHORS_TO_REFERENCE {
    label 'liftoff'
    tag "${sample}:${chromosome}"
    publishDir "${params.outdir}/06_lift_to_ref/liftoff", mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(chromosome), path(query_fasta), path(anchor_gff),
          path(anchor_map), path(features), path(reference_fasta)

    output:
    tuple val(sample), val(chromosome),
          path("${sample}_${chromosome}_liffout.gff_polished"), emit: polished
    tuple val(sample), val(chromosome),
          path("${sample}_${chromosome}.liftoff.log"), emit: logs

    script:
    """
    set -euo pipefail
    cp -L '${reference_fasta}' target.fa
    cp -L '${query_fasta}' reference.fa
    cp '${anchor_gff}' annotation.gff
    cp '${features}' features.txt
    mkdir -p liftoff_tmp
    if ! grep -qvE '^[[:space:]]*(#|\$)' annotation.gff; then
        printf 'INFO: no anchors for ${sample}:${chromosome}; Liftoff was not required.\n' \
            > '${sample}_${chromosome}.liftoff.log'
        printf '##gff-version 3\n' > '${sample}_${chromosome}_liffout.gff_polished'
    else
        set +e
        liftoff -g annotation.gff \
            -o '${sample}_${chromosome}_liffout.gff' \
            -exclude_partial -dir liftoff_tmp -polish -copies \
            -f features.txt -overlap 1 -p ${task.cpus} \
            target.fa reference.fa \
            > '${sample}_${chromosome}.liftoff.log' 2>&1
        status=\$?
        set -e
        if [[ \${status} -ne 0 ]]; then
            echo "ERROR: anchor Liftoff failed for ${sample}:${chromosome} (status \${status})" >&2
            tail -n 50 '${sample}_${chromosome}.liftoff.log' >&2 || true
            exit \${status}
        fi
        if [[ ! -s '${sample}_${chromosome}_liffout.gff_polished' ]]; then
            echo "ERROR: anchor Liftoff finished without a polished GFF for ${sample}:${chromosome}" >&2
            tail -n 50 '${sample}_${chromosome}.liftoff.log' >&2 || true
            exit 1
        fi
        python '${projectDir}/bin/burrite_stage.py' restore-anchor-ids \
            --gff '${sample}_${chromosome}_liffout.gff_polished' \
            --mapping '${anchor_map}' \
            --output '${sample}_${chromosome}_liffout.gff_polished.restored'
        mv '${sample}_${chromosome}_liffout.gff_polished.restored' \
            '${sample}_${chromosome}_liffout.gff_polished'
    fi
    """
}


process BUILD_SAMPLE_BED {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/06_lift_to_ref", mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.endsWith('.bed') ? "bed_files/${filename}" : filename }

    input:
    tuple val(sample), path(polished_gffs)

    output:
    tuple val(sample), path("${sample}.bed"), path("transfer_${sample}"), emit: bed

    script:
    def gffOptions = polished_gffs.collect { "--gff '${it}'" }.join(' ')
    """
    python '${projectDir}/bin/burrite_stage.py' sample-bed \
        --sample '${sample}' \
        ${gffOptions}
    """
}


process REFERENCE_BED {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/06_lift_to_ref/bed_files", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(annotation)
    path reference

    output:
    tuple val(sample), path("${sample}.bed"), emit: bed

    script:
    def maxOption = params.max_length != null ? "--max-length ${params.max_length}" : ''
    """
    python '${projectDir}/bin/burrite_stage.py' reference-bed \
        --sample '${sample}' \
        --reference '${reference}' \
        --annotation '${annotation}' \
        --min-length ${params.min_length} \
        ${maxOption}
    """
}


process VALIDATE_SAMPLE_BEDS {
    label 'small'
    tag 'per-genome BED validation'
    publishDir "${params.outdir}/06_lift_to_ref/validation", mode: 'copy', overwrite: true

    input:
    path beds
    val expected_samples
    val allow_empty

    output:
    path 'validated_beds/*.bed', emit: beds
    path 'sample_bed_validation.tsv', emit: report

    script:
    def bedList = beds instanceof List ? beds : [beds]
    def bedOptions = bedList.collect { "--bed '${it}'" }.join(' ')
    def allowOption = enabled(allow_empty) ? '--allow-empty' : ''
    """
    python '${projectDir}/bin/burrite_stage.py' validate-sample-beds \
        ${bedOptions} \
        --expected '${expected_samples}' \
        --output-dir validated_beds \
        --report sample_bed_validation.tsv \
        ${allowOption}
    """
}


process COLLAPSE_LOCI {
    label 'small'
    tag "${params.merge_method}:${params.merge_by}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path beds

    output:
    path '07_merge', emit: merge_dir

    script:
    def bedOptions = beds.collect { "--bed '${it}'" }.join(' ')
    """
    export PYTHONHASHSEED=0
    python '${projectDir}/bin/burrite_stage.py' collapse \
        ${bedOptions} \
        --method '${params.merge_method}' \
        --merge-by '${params.merge_by}' \
        --lenthr ${params.lenthr} \
        --minlen ${params.minlen} \
        --bedmin ${params.bedmin}
    """
}


process GENERATE_CONSENSUS_FLANKS {
    label 'small'
    tag 'consensus flanks'
    publishDir "${params.outdir}/07_merge", mode: 'copy', overwrite: true

    input:
    path merge_dir
    path chromosome_list

    output:
    path 'FM_*.gff', emit: gffs
    path '*_gff_features.txt', emit: features

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' consensus-flanks \
        --collapsed-bed '${merge_dir}/collapsed.bed' \
        --chromosome-list '${chromosome_list}' \
        --flank-size ${params.flank_size}
    """
}


process LIFT_CONSENSUS_TO_ASSEMBLY {
    label 'liftoff'
    tag "${sample}:${chromosome}"
    publishDir "${params.outdir}/07_merge/liftback_liftoff", mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(chromosome), path(query_fasta), path(consensus_gff),
          path(features), path(reference_fasta)

    output:
    tuple val(sample), val(chromosome),
          path("${sample}_${chromosome}_consensus_liffout.gff_polished"), emit: polished
    tuple val(sample), val(chromosome),
          path("${sample}_${chromosome}.consensus_liftoff.log"), emit: logs

    script:
    """
    set -euo pipefail
    cp -L '${query_fasta}' target.fa
    cp -L '${reference_fasta}' reference.fa
    cp '${consensus_gff}' annotation.gff
    cp '${features}' features.txt
    mkdir -p liftoff_tmp
    if ! grep -qvE '^[[:space:]]*(#|\$)' annotation.gff; then
        printf 'INFO: no consensus anchors for ${sample}:${chromosome}; Liftoff was not required.\n' \
            > '${sample}_${chromosome}.consensus_liftoff.log'
        printf '##gff-version 3\n' > '${sample}_${chromosome}_consensus_liffout.gff_polished'
    else
        set +e
        liftoff -g annotation.gff \
            -o '${sample}_${chromosome}_consensus_liffout.gff' \
            -exclude_partial -dir liftoff_tmp -polish -copies \
            -f features.txt -overlap 1 -p ${task.cpus} \
            target.fa reference.fa \
            > '${sample}_${chromosome}.consensus_liftoff.log' 2>&1
        status=\$?
        set -e
        if [[ \${status} -ne 0 ]]; then
            echo "ERROR: consensus Liftoff failed for ${sample}:${chromosome} (status \${status})" >&2
            tail -n 50 '${sample}_${chromosome}.consensus_liftoff.log' >&2 || true
            exit \${status}
        fi
        if [[ ! -s '${sample}_${chromosome}_consensus_liffout.gff_polished' ]]; then
            echo "ERROR: consensus Liftoff finished without a polished GFF for ${sample}:${chromosome}" >&2
            tail -n 50 '${sample}_${chromosome}.consensus_liftoff.log' >&2 || true
            exit 1
        fi
    fi
    """
}


process EVALUATE_SAMPLE {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/08_liftback", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(polished_gffs)
    path merge_dir

    output:
    tuple val(sample),
          path("${sample}_Transfers_TE.txt"),
          path("${sample}_liffout_polished.gff"), emit: transfers

    script:
    def gffOptions = polished_gffs.collect { "--gff '${it}'" }.join(' ')
    """
    python '${projectDir}/bin/burrite_stage.py' evaluate-sample \
        --sample '${sample}' \
        --collapsed-bed '${merge_dir}/collapsed.bed' \
        ${gffOptions}
    """
}


process EVALUATE_REFERENCE {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/08_liftback", mode: 'copy', overwrite: true

    input:
    val sample
    path merge_dir

    output:
    tuple val(sample), path("${sample}_Transfers_TE.txt"), emit: transfer

    script:
    """
    python '${projectDir}/bin/burrite_stage.py' evaluate-reference \
        --sample '${sample}' \
        --collapsed-bed '${merge_dir}/collapsed.bed' \
        --flank-size ${params.flank_size}
    """
}


process VALIDATE_TRANSFER_FILES {
    label 'small'
    tag 'per-genome transfer validation'
    publishDir "${params.outdir}/08_liftback/validation", mode: 'copy', overwrite: true

    input:
    path transfer_files
    path merge_dir
    val expected_samples

    output:
    path 'validated_transfers/*_Transfers_TE.txt', emit: transfers
    path 'transfer_validation.tsv', emit: report

    script:
    def transferList = transfer_files instanceof List ? transfer_files : [transfer_files]
    def transferOptions = transferList.collect { "--transfer '${it}'" }.join(' ')
    """
    python '${projectDir}/bin/burrite_stage.py' validate-transfers \
        ${transferOptions} \
        --expected '${expected_samples}' \
        --collapsed-bed '${merge_dir}/collapsed.bed' \
        --output-dir validated_transfers \
        --report transfer_validation.tsv
    """
}


process BUILD_MATRIX {
    label 'small'
    tag 'presence/absence matrix'
    publishDir "${params.outdir}/09_final", mode: 'copy', overwrite: true

    input:
    path merge_dir
    path transfer_files
    val sample_order

    output:
    path 'Presence_Absence_matrix.txt', emit: matrix

    script:
    def transferOptions = transfer_files.collect { "--transfer '${it}'" }.join(' ')
    """
    export PYTHONHASHSEED=0
    python '${projectDir}/bin/burrite_stage.py' matrix \
        --collapsed-bed '${merge_dir}/collapsed.bed' \
        --merge-dir '${merge_dir}' \
        ${transferOptions} \
        --samples '${sample_order}' \
        --method '${params.merge_method}'
    """
}


process REFERENCE_RENAME_SENTINEL {
    label 'small'
    tag "${sample}"

    input:
    val sample

    output:
    tuple val(sample), path("${sample}.identity.fasta"), emit: renamed

    script:
    """
    : > '${sample}.identity.fasta'
    """
}


process FINAL_GENOME_GFF {
    label 'small'
    tag "${sample}"
    publishDir "${params.outdir}/09_final/per_genome_gff", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(original_fasta), path(annotation), path(filtered_bed), path(renamed_fasta)
    path matrix
    path merge_dir
    path transfer_files

    output:
    tuple val(sample), path("${sample}*TEs.gff"), emit: gff

    script:
    def transferOptions = transfer_files.collect { "--transfer '${it}'" }.join(' ')
    def maxOption = params.max_length != null ? "--max-length ${params.max_length}" : ''
    """
    export PYTHONHASHSEED=0
    python '${projectDir}/bin/burrite_stage.py' final-gff \
        --sample '${sample}' \
        --matrix '${matrix}' \
        --merge-dir '${merge_dir}' \
        ${transferOptions} \
        --original-fasta '${original_fasta}' \
        --renamed-fasta '${renamed_fasta}' \
        --annotation '${annotation}' \
        --filtered-bed '${filtered_bed}' \
        --mode '${params.final_annotation}' \
        --min-length ${params.min_length} \
        ${maxOption} \
        --overlap ${params.final_annotation_dedup_overlap}
    """
}


workflow {
    if( !params.assemblies || !params.te_lib || !params.reference || !params.chromosome_list ) {
        error "Required parameters: --assemblies --te_lib --reference --chromosome_list"
    }
    if( !['incremental', 'groupby'].contains(params.merge_method.toString()) ) {
        error "--merge_method must be incremental or groupby"
    }
    if( !['family', 'superfamily'].contains(params.merge_by.toString()) ) {
        error "--merge_by must be family or superfamily"
    }
    if( !['burrite', 'all'].contains(params.final_annotation.toString()) ) {
        error "--final_annotation must be burrite or all"
    }
    if( !['Y', 'YES', 'N', 'NO', 'TRUE', 'FALSE', '1', '0'].contains(params.include_reference.toString().toUpperCase()) ) {
        error "--include_reference must be Y or N"
    }
    if( (params.flank_size as int) < 1 || (params.min_length as int) < 0 ||
        (params.lenthr as int) < 0 || (params.minlen as int) < 0 || (params.bedmin as int) < 0 ) {
        error "Length parameters must be non-negative and --flank_size must be at least 1"
    }
    if( params.max_length != null && (params.max_length as int) < (params.min_length as int) ) {
        error "--max_length must be greater than or equal to --min_length"
    }
    if( (params.final_annotation_dedup_overlap as double) < 0.0 ||
        (params.final_annotation_dedup_overlap as double) > 1.0 ) {
        error "--final_annotation_dedup_overlap must be between 0 and 1"
    }

    def assembliesFile = file(params.assemblies, checkIfExists: true)
    def teLibrary = file(params.te_lib, checkIfExists: true)
    def reference = file(params.reference, checkIfExists: true)
    def chromosomeList = file(params.chromosome_list, checkIfExists: true)
    def onecodeDirectory = file(params.onecode_dir, checkIfExists: true)
    // v3 deliberately uses the token before the first dot so a reference such
    // as 100002.Chr_scaffolds.fa is called 100002, like GraffiTE's sample IDs.
    def referenceName = reference.getFileName().toString().tokenize('.')[0]
    def includeReference = enabled(params.include_reference)

    def onecodeScripts = [
        onecodeDirectory.resolve('build_dictionary.pl'),
        onecodeDirectory.resolve('one_code_to_find_them_all_but_sanely.pl')
    ]
    def missingOnecode = onecodeScripts.findAll { !java.nio.file.Files.isRegularFile(it) }
    if( missingOnecode ) {
        error "Missing BurriTE-compatible OneCode scripts: ${missingOnecode.join(', ')}"
    }

    def rows = java.nio.file.Files.readAllLines(assembliesFile).drop(1).findAll { it.trim() }
    def alternativeRecords = []
    rows.eachWithIndex { line, index ->
        def columns = line.split(',', -1)
        if( columns.size() < 2 ) {
            error "Invalid assemblies CSV row ${index + 2}: expected FASTA,sample"
        }
        def sample = columns[1].trim()
        if( !(sample ==~ /[A-Za-z0-9_-]+/) || sample.contains('.') ) {
            error "Invalid sample ID '${sample}' at assemblies CSV row ${index + 2}"
        }
        alternativeRecords << tuple(sample, resolveCsvPath(columns[0].trim(), assembliesFile.getParent()))
    }
    def alternativeNames = alternativeRecords.collect { it[0] }
    if( !alternativeRecords ) {
        error 'The assemblies CSV must contain at least one alternative assembly'
    }
    if( alternativeNames.size() != alternativeNames.toSet().size() ) {
        error 'The assemblies CSV contains duplicated sample IDs'
    }
    if( alternativeNames.contains(referenceName) ) {
        error "The reference ID '${referenceName}' is also present in the alternative assemblies CSV"
    }

    def chromosomes = java.nio.file.Files.readAllLines(chromosomeList)
        .collect { it.trim() }.findAll { it }
    if( !chromosomes ) {
        error 'The chromosome list is empty'
    }

    def allRecords = alternativeRecords + [tuple(referenceName, reference)]
    def outputRoot = java.nio.file.Paths.get(params.outdir.toString())
    if( !outputRoot.isAbsolute() ) {
        outputRoot = java.nio.file.Paths.get(launchDir.toString()).resolve(outputRoot)
    }
    outputRoot = outputRoot.normalize().toAbsolutePath()
    def graffiteTmp = params.graffite_tmpdir \
        ? java.nio.file.Paths.get(params.graffite_tmpdir.toString()).toAbsolutePath().normalize() \
        : outputRoot.resolve('01_GraffiTE/tmp')
    def graffiteImage = params.graffite_image \
        ? file(params.graffite_image, checkIfExists: true).toAbsolutePath().toString() : ''

    alt_assemblies_ch = Channel.fromList(alternativeRecords)
    all_assemblies_ch = Channel.fromList(allRecords)
    te_library_ch = Channel.value(teLibrary)
    reference_ch = Channel.value(reference)
    chromosome_list_ch = Channel.value(chromosomeList)
    onecode_dir_ch = Channel.value(onecodeDirectory)

    if( params.graffite_vcf ) {
        graffite_vcf_ch = Channel.value(file(params.graffite_vcf, checkIfExists: true))
    }
    else {
        RUN_GRAFFITE(
            Channel.value(assembliesFile), te_library_ch, reference_ch,
            Channel.value(outputRoot.toString()), Channel.value(graffiteTmp.toString()),
            Channel.value(graffiteImage)
        )
        graffite_vcf_ch = RUN_GRAFFITE.out.vcf
    }

    PARSE_GRAFFITE_VCF(
        graffite_vcf_ch,
        Channel.value(alternativeNames.join(',')),
        Channel.value(referenceName)
    )
    graffite_rm_ch = PARSE_GRAFFITE_VCF.out.records
        .flatten()
        .map { record -> tuple(record.baseName.replaceFirst(/_GraffiTE$/, ''), record) }

    REPEATMASKER(all_assemblies_ch, te_library_ch)
    ONECODE(REPEATMASKER.out.rm_dir, onecode_dir_ch)
    onecode_annotation_ch = ONECODE.out.result.map { sample, annotation, dictionary, log ->
        tuple(sample, annotation)
    }

    conversion_input_ch = onecode_annotation_ch.join(graffite_rm_ch)
    CONVERT_ANNOTATIONS(conversion_input_ch)
    rm_gff_ch = CONVERT_ANNOTATIONS.out.gffs.map { sample, rm, graffite -> tuple(sample, rm) }
    graffite_gff_ch = CONVERT_ANNOTATIONS.out.gffs.map { sample, rm, graffite -> tuple(sample, graffite) }

    alt_graffite_input_ch = Channel.fromList(alternativeRecords)
        .join(graffite_gff_ch.filter { sample, gff -> sample != referenceName })
    LIFT_GRAFFITE(alt_graffite_input_ch, reference_ch)
    reference_graffite_input_ch = graffite_gff_ch.filter { sample, gff -> sample == referenceName }
    COPY_REFERENCE_GRAFFITE(reference_graffite_input_ch)
    lifted_graffite_ch = LIFT_GRAFFITE.out.lifted.mix(COPY_REFERENCE_GRAFFITE.out.lifted)

    JOIN_ANNOTATIONS(rm_gff_ch.join(lifted_graffite_ch))
    FILTER_ANNOTATION(JOIN_ANNOTATIONS.out.annotation)
    MAKE_ANCHORS(FILTER_ANNOTATION.out.bed)

    MINIMAP_ASSEMBLY(Channel.fromList(alternativeRecords), reference_ch)
    rename_input_ch = Channel.fromList(alternativeRecords)
        .join(MAKE_ANCHORS.out.anchors)
        .join(MINIMAP_ASSEMBLY.out.paf)
    RENAME_ASSEMBLY(rename_input_ch)

    EXTRACT_REFERENCE_CHROM(Channel.fromList(chromosomes), reference_ch)
    alt_chromosome_input_ch = RENAME_ASSEMBLY.out.renamed
        .combine(Channel.fromList(chromosomes))
        .map { sample, renamed, anchors, anchor_map, chromosome ->
            tuple(sample, chromosome, renamed, anchors, anchor_map)
        }
    PREPARE_ALT_CHROM(alt_chromosome_input_ch)

    anchor_jobs_ch = PREPARE_ALT_CHROM.out.prepared
        .map { sample, chromosome, query, anchors, anchor_map, features ->
            tuple(chromosome, sample, query, anchors, anchor_map, features)
        }
        // Multiple assemblies have the same chromosome key. `combine(by: 0)`
        // deliberately creates one job per assembly/chromosome pair; `join`
        // would retain only one duplicate key and silently drop assemblies.
        .combine(EXTRACT_REFERENCE_CHROM.out.chromosome_fasta, by: 0)
        .map { chromosome, sample, query, anchors, anchor_map, features, refchr ->
            tuple(sample, chromosome, query, anchors, anchor_map, features, refchr)
        }
    LIFT_ANCHORS_TO_REFERENCE(anchor_jobs_ch)
    sample_anchor_groups_ch = LIFT_ANCHORS_TO_REFERENCE.out.polished
        .map { sample, chromosome, gff -> tuple(sample, gff) }
        .groupTuple()
    BUILD_SAMPLE_BED(sample_anchor_groups_ch)

    alternative_beds_ch = BUILD_SAMPLE_BED.out.bed.map { sample, bed, intermediates -> tuple(sample, bed) }
    if( includeReference ) {
        reference_annotation_ch = JOIN_ANNOTATIONS.out.annotation
            .filter { sample, annotation -> sample == referenceName }
        REFERENCE_BED(reference_annotation_ch, reference_ch)
        collapse_beds_ch = alternative_beds_ch.mix(REFERENCE_BED.out.bed)
        bedNames = [referenceName] + alternativeNames
    }
    else {
        collapse_beds_ch = alternative_beds_ch
        bedNames = alternativeNames
    }

    VALIDATE_SAMPLE_BEDS(
        collapse_beds_ch.map { sample, bed -> bed }.collect(),
        Channel.value(bedNames.join(',')),
        Channel.value(params.allow_empty_sample_beds)
    )
    COLLAPSE_LOCI(VALIDATE_SAMPLE_BEDS.out.beds)
    merge_dir_value = COLLAPSE_LOCI.out.merge_dir
    GENERATE_CONSENSUS_FLANKS(merge_dir_value, chromosome_list_ch)

    consensus_by_chromosome_ch = GENERATE_CONSENSUS_FLANKS.out.gffs
        .flatten()
        .map { gff -> tuple(gff.baseName.replaceFirst(/^FM_/, ''), gff) }
        .join(
            GENERATE_CONSENSUS_FLANKS.out.features.flatten().map { features ->
                tuple(features.baseName.replaceFirst(/_gff_features$/, ''), features)
            }
        )
        .join(EXTRACT_REFERENCE_CHROM.out.chromosome_fasta)

    consensus_jobs_ch = PREPARE_ALT_CHROM.out.prepared
        .map { sample, chromosome, query, anchors, anchor_map, features ->
            tuple(chromosome, sample, query)
        }
        // As above, expand each chromosome consensus to every assembly.
        .combine(consensus_by_chromosome_ch, by: 0)
        .map { chromosome, sample, query, consensus, features, refchr ->
            tuple(sample, chromosome, query, consensus, features, refchr)
        }
    LIFT_CONSENSUS_TO_ASSEMBLY(consensus_jobs_ch)
    consensus_groups_ch = LIFT_CONSENSUS_TO_ASSEMBLY.out.polished
        .map { sample, chromosome, gff -> tuple(sample, gff) }
        .groupTuple()
    EVALUATE_SAMPLE(consensus_groups_ch, merge_dir_value)

    alternative_transfer_ch = EVALUATE_SAMPLE.out.transfers.map { sample, table, combined ->
        tuple(sample, table)
    }
    if( includeReference ) {
        EVALUATE_REFERENCE(Channel.value(referenceName), merge_dir_value)
        all_transfer_ch = alternative_transfer_ch.mix(EVALUATE_REFERENCE.out.transfer)
        matrixNames = [referenceName] + alternativeNames
    }
    else {
        all_transfer_ch = alternative_transfer_ch
        matrixNames = alternativeNames
    }
    raw_transfer_files_value = all_transfer_ch.map { sample, table -> table }.collect()
    VALIDATE_TRANSFER_FILES(
        raw_transfer_files_value,
        merge_dir_value,
        Channel.value(matrixNames.join(','))
    )
    transfer_files_value = VALIDATE_TRANSFER_FILES.out.transfers
    BUILD_MATRIX(merge_dir_value, transfer_files_value, Channel.value(matrixNames.join(',')))
    matrix_value = BUILD_MATRIX.out.matrix

    renamed_alt_ch = RENAME_ASSEMBLY.out.renamed.map { sample, renamed, anchors, anchor_map ->
        tuple(sample, renamed)
    }
    if( includeReference ) {
        REFERENCE_RENAME_SENTINEL(Channel.value(referenceName))
        renamed_for_final_ch = renamed_alt_ch.mix(REFERENCE_RENAME_SENTINEL.out.renamed)
        finalSourceRecords = allRecords
    }
    else {
        renamed_for_final_ch = renamed_alt_ch
        finalSourceRecords = alternativeRecords
    }

    final_metadata_ch = Channel.fromList(finalSourceRecords)
        .join(JOIN_ANNOTATIONS.out.annotation)
        .join(FILTER_ANNOTATION.out.bed)
        .join(renamed_for_final_ch)
    FINAL_GENOME_GFF(final_metadata_ch, matrix_value, merge_dir_value, transfer_files_value)
}
