from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from cov2019.utils import srr_to_url
HTTP = HTTPRemoteProvider()

genomes=config['genomes']
samples_file=config['samples']
threads=9
samples=list()

with open(samples_file) as fin:
    samples = [sample.rstrip('\r\n') for sample in fin]

#print(samples)

def get_reads(wildcards):
    link=srr_to_url(wildcards.sample)
    return HTTP.remote(link)

rule all:
    input:
        #expand("{sample}.fasta", sample=samples)
        #expand("{sample}.1_{pair}.fastq", sample=samples, pair=[1,2])
        expand("montagens/{sample}.fa", sample=samples)


rule fastq_dump:
    input:
        reads=get_reads
    params:
        sra="data/{sample}"
    output:
        "data/{sample}_1.fastq",
        "data/{sample}_2.fastq"
    conda:
        "envs/filter.yaml"
    log:
        "fastq-dump_{sample}.log"
    shell:
        "mv {input.reads} {params.sra};fastq-dump -O data --split-files {params.sra}"

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
        "bwa index -p {output} {input.genomes}; touch ref/index"

rule filter_virus_reads:
    input:
        read1="data/{sample}_1.fastq",
        read2="data/{sample}_2.fastq",
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
        "samtools view -bS -f 2 - | "
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
        "spades.py -t {threads} -1 {input[0]} -2 {input[1]} --careful  -o {output}"
