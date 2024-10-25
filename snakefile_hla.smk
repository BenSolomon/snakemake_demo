# snakemake -s snakefile_hla.smk --use-conda
# snakemake -s snakefile_hla.smk --dag | dot -Tpng > hla_dag.png

rule all:
    input:
        # expand('hisat/sample_{num}.bam', num=[1,2]), # HISAT output, note this is not needed when downstream outputs are requested
        expand('arcasHLA/sample_{num}.genotype.json', num=[1,2]), # arcasHLA output
        expand('optitype/sample_{num}_result.tsv', num=[1,2]), # optitype output
        expand('phlat/sample_{num}_HLA.sum', num=[1,2]) # phlat output

rule hisat:
    input:
        "{sample}.fq.gz"
    output:
        "hisat/{sample}.bam",
        "hisat/{sample}.bam.bai"
    conda:
        "samtools"
    threads: 2
    params:
        GENOME_INDEX='/labs/khatrilab/solomonb/references/grch38/hisat/genome'
    shell:
        '''
        hisat2 --dta -x {params.GENOME_INDEX} \
        -U {input} |
        samtools sort -@ {threads} -O BAM > {output[0]}

        samtools index -@ {threads} {output[0]} {output[1]}
        '''

rule arcasHLA:
    input:
        "hisat/{sample}.bam"
    output:
        "arcasHLA/{sample}.extracted.fq.gz",
        "arcasHLA/{sample}.genes.json",
        "arcasHLA/{sample}.genotype.json",
        temp("arcasHLA/{sample}.genotype.log"),
        temp("arcasHLA/{sample}.alignment.p"),
        temp("arcasHLA/{sample}.extract.log")
    conda:
        "samtools"
    threads: 2
    params:
        ARCAS_BIN='/labs/khatrilab/solomonb/software_server/arcasHLA/arcasHLA',
        ARCAS_DIR='arcasHLA/'
    shell:
        '''
        {params.ARCAS_BIN} extract {input} --unmapped -t {threads} -o {params.ARCAS_DIR}

        {params.ARCAS_BIN} genotype {output[0]} -t {threads} -o {params.ARCAS_DIR}
        '''

rule optitype:
    input:
        "arcasHLA/{sample}.extracted.fq.gz"
    output:
        "optitype/{sample}_result.tsv",
        "optitype/{sample}_coverage_plot.pdf"
    conda:
        "optitype"
    params:
        OPTITYPE_DIR='optitype/'
    shell:
        '''
        OptiTypePipeline.py -i {input} \
        -o {params.OPTITYPE_DIR} \
        -p {wildcards.sample} \
        --rna \
        -v
        '''

rule phlat:
    input:
        "arcasHLA/{sample}.extracted.fq.gz"
    output:
        "phlat/{sample}_HLA.sum"
    conda:
        "phlat"
    threads: 2
    params:
        PHLAT_BIN_DIR='/labs/khatrilab/solomonb/software_server/phlat-release',
        PHLAT_BIN='/labs/khatrilab/solomonb/software_server/phlat-release/dist/PHLAT.py',
        PHLAT_INDEX='/labs/khatrilab/solomonb/software_server/phlat-release/b2folder',
        BOWTIE_BIN='/labs/khatrilab/solomonb/software_server/bowtie2/bowtie2-2.0.0-beta7/bowtie2',
        PHLAT_DIR='phlat'
    shell:
        '''
        python -O {params.PHLAT_BIN} \
        -1 {input} \
        -index {params.PHLAT_INDEX} \
        -b2url {params.BOWTIE_BIN} \
        -tag {wildcards.sample} \
        -e {params.PHLAT_BIN_DIR} \
        -o {params.PHLAT_DIR} \
        -p {threads} \
        -pe 0 \
        -tmp 0 
        '''


