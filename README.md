# cov2019
Cross-species RNA analysis between BetaCoV 2019-2020 and animals

Sequencia de coronavirus (GenomicFastaResults.fasta):
https://www.viprbrc.org/brc/home.spg?decorator=corona

Linha de comando do cd-hit:

`cd-hit-est -i GenomicFastaResults.fasta -o GenomicFastaResults.cdhit98.fasta -g 1 -c 0.98 -M 4000 -n 6`
