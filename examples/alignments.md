# seals-2 working with alignments
## Preparing an alignment for converting to HMM or PSSM
### Source data:
### Source data:
`protset.afa` aligned protein sequences in FASTA format.
### Command:
```
fa2sr -w=0 protset.afa | sr_filter -conplus -hcon=0 -gcon=0.5 -ncon='protset' -grcut=0.5 | sr2fa > profile.FASTA
```
converts the `protset.afa` alignment into the "seqrows" format (`prot_id sequence`), adds consensus to sites [with](https://doi.org/10.1093/ve/veab015) the (weighted) gap content no more than 0.5 and name "protset"; removes sites with the (weighted) gap content over 0.5, converts the alignment back to FASTA format and saves `profile.FASTA` file.

## Preparing an alignment for tree reconstruction
### Source data:
`protset.afa` aligned protein sequences in FASTA format.

`prot_genome.tab` of the following tab-delimited format:
```
prot_id genome_id
```
### Commands:
```
fa2sr -w=0 protset.afa | sr_scoreMatch | tab_select -k1=1 -kcut=0.2 | tab_merge -t= prot_genome.tab | tab_best -k1=3 -k2=1 -m=2 > tmp.tab
```
converts the `protset.afa` alignment into the "seqrows" format (`prot_id sequence`), assigns pairwise BLOSUM62 scores to each sequence against the consensus, discards lines with relative scores below 0.2, adds genome information to the table, selects the highest-scoring sequence for each genome and prints out the resulting score table.
```
fa2sr -w=0 protset.afa | tab_select -t= tmp.tab -k1=1 -k2=2 | sr_filter -grcut=0.5 -hocut= 0.1 | sr2fa > prot_for_tree.afa
```
converts the `protset.afa` alignment into the "seqrows" format (`prot_id sequence`), selects sequences matching the `tmp.tab` table, removes sites [with](https://doi.org/10.1093/ve/veab015) the (weighted) gap content over 0.5 and homogeneity below 0.1, converts the alignment back to FASTA format and saves `prot_for_tree.afa` file.
