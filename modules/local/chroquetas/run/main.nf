// ChroQueTaS — antifungal-resistance scan (FungAMR panel). Mirrors steps/15_chroquetas.sh.
// Channel is nmquijada (not bioconda); no nf-core module exists.
process CHROQUETAS {
    tag "${meta.id}"
    label 'process_medium'

    conda "nmquijada::chroquetas=1.0.1"

    input:
    tuple val(meta), path(assembly), val(species_arg)

    output:
    tuple val(meta), path("${meta.id}.ChroQueTaS.AMR_summary.txt"), emit: summary, optional: true
    tuple val(meta), path("${meta.id}.chroquetas.log"),             emit: log
    tuple val(meta), path("SKIPPED.txt"),                           emit: skipped, optional: true
    tuple val("${task.process}"), val('chroquetas'), eval('ChroQueTas.sh --version 2>&1 | head -1 | sed -r "s/\\x1b\\[[0-9;]*[a-zA-Z]//g; s/\\x1b\\([A-Z]//g; s/.*v//"'), topic: versions, emit: versions_chroquetas

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    if [ -z "${species_arg}" ] || [ "${species_arg}" = "NA" ]; then
        echo "skipped: species not supported by ChroQueTas 1.0.1" > SKIPPED.txt
        touch ${meta.id}.chroquetas.log
        exit 0
    fi
    ChroQueTas.sh -g ${assembly} -o ${meta.id}_chroquetas -s ${species_arg} \\
        -t ${task.cpus} --trans_code 12 2> ${meta.id}.chroquetas.log
    cp ${meta.id}_chroquetas/${meta.id}.ChroQueTaS.AMR_summary.txt .
    """

    stub:
    """
    printf "Protein\\tFragment\\tPosition_reference\\tAA_reference\\tAA_query\\tFungicide_resistance\\n" > ${meta.id}.ChroQueTaS.AMR_summary.txt
    touch ${meta.id}.chroquetas.log
    """
}
