# seals-2 working with tables
## Taxonomy of best BLAST hits
### Source data:
`blasthits.tab` of the following tab-delimited format:
```
query_id subject_id ... score
```
in columns \## 1, 2 and 12 (`-outfmt 6` in BLAST command line).

`prot_tax.tab` of the following tab-delimited format:
```
subject taxon
```
in columns \## 1 and 2.
### Command:
```
tab_best blasthits.tab -k1=1 -k2=12 -m=2 | tab_shuffle -l='1,2' | tab_merge -t= prot_tax.tab -k1=2 -k2=1 | tab_aggr -k1=3 | sort -k 2,2nr -k 1,1
```
will select the highest-scoring hit for each query from `blasthits.tab`; trim the line to `query_id subject_id` pair; add taxon data to every line, then count the number of times each taxon is present among the best hits and print out the sorted list of taxa
