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
        expand("{sample}.fa", sample='SRR1869462')


rule fastq_dump:
    input:
        reads=get_reads
    output:
        "{sample}.1_1.fastq",
        "{sample}.1_2.fastq"
    conda:
        "envs/filter.yaml"
    log:

    shell:
        "fastq-dump --split-files {input.reads}"

rule build_index:
    input:
        genomes=genomes
    output:
        "index"
    conda:
        "envs/filter.yaml"
    log:
        "logs/index.log"
    shell:
        "touch index;bwa index -p {output} {input.genomes}"

rule filter_virus_reads:
    input:
        read1="{sample}.1_1.fastq",
        read2="{sample}.1_2.fastq",
        index="index"
    output:
        "{sample}_1.fq",
        "{sample}_2.fq"
    conda:
        "envs/filter.yaml"
    log:
        "log/filter_{sample}.log"
    threads: threads
    shell:
        "bwa mem -t {threads} {input.index} {input.read1} {input.read2} | "
        "samtools view -bS - | "
        "bamToFastq -i /dev/stdin -fq {output[0]} -fq2 {output[1]}"

rule assembly_virus:
    input:
        "{sample}_1.fq",
        "{sample}_2.fq"
    output:
        "{sample}.fa"
    conda:
        "envs/assembly.yaml"
    log:
        "logs/assembly_{sample}.log"
    threads:threads
    shell:
        "spades.py -1 {input[0]} -2 {input[1]} --careful  -o {output}"

