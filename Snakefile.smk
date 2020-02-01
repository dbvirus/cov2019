from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

genomes=config['genomes']
samples_file=config['samples']
threads=10
samples=list()

with open(samples_file) as fin:
    samples = [sample.rstrip('\r\n') for sample in fin]

#print(samples)

def get_reads(wildcards):
    
    return HTTP.remote("https://sra-download.st-va.ncbi.nlm.nih.gov/sos2/sra-pub-run-6/SRR1869462/SRR1869462.1", keep_local=True)

rule all:
    input:
        #expand("{sample}.fasta", sample=samples)
        #expand("{sample}.1_{pair}.fastq", sample=samples, pair=[1,2])
        expand("montagens/{sample}.fa", sample='SRR1869462')


rule fastq_dump:
    input:
        reads=get_reads
    output:
        "data/{sample}.1_1.fastq",
        "data/{sample}.1_2.fastq"
    conda:
        "envs/filter.yaml"
    log:

    shell:
        "cd data/;fastq-dump --split-files {input.reads}"

rule build_index:
    input:
        genomes=genomes
    output:
        "ref/index"
    conda:
        "envs/filter.yaml"
    log:
        "logs/index.log"
    shell:
        "mkdir ref;touch ref/index;bwa index -p {output} {input.genomes}"

rule filter_virus_reads:
    input:
        read1="data/{sample}.1_1.fastq",
        read2="data/{sample}.1_2.fastq",
        index="ref/index"
    output:
        "data/{sample}_1.fq",
        "data/{sample}_2.fq"
    conda:
        "envs/filter.yaml"
    log:
        "log/filter_{sample}.log"
    threads: threads
    shell:
        "bwa mem -t {threads} {input.index} {input.read1} {input.read2} | "
        "samtools view -bS -F 4 - | "
        "bamToFastq -i /dev/stdin -fq {output[0]} -fq2 {output[1]}"

rule assembly_virus:
    input:
        "data/{sample}_1.fq",
        "data/{sample}_2.fq"
    output:
        "montagens/{sample}.fa"
    conda:
        "envs/assembly.yaml"
    log:
        "logs/assembly_{sample}.log"
    threads:threads
    shell:
        "spades.py -1 {input[0]} -2 {input[1]} --careful  -o {output}"

