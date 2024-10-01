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
subject_id taxon
```
in columns \## 1 and 2.
### Command:
```
tab_best blasthits.tab -k1=1 -k2=12 -m=2 | tab_shuffle -l='1,2' | tab_merge -t= prot_tax.tab -k1=2 -k2=1 | tab_aggr -k1=3 | sort -k 2,2nr -k 1,1
```
selects the highest-scoring hit for each query from `blasthits.tab`; trims the line to `query_id subject_id` pair; adds taxon data to every line, then counts the number of times each taxon is present among the best hits and prints out the sorted list of taxa counts.

## Annotate a chromosome with COGs
### Source data (https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/ and GenBank):
`cog-24.cog.csv` of the following comma-delimited format:
```
gene_id,...,...,...,cog_loc,...,cog_id,...
```
in columns \## 1, 5 and 7.

`chromosome.lst` containing a list of `gene_id` in the order, encoded on the chromosome.
### Commands:
```
tab_select cog-24.cog.csv -t= chromosome.lst -k1=1 -s1=',' | sort -k 1,1 -k 5,5n -t ',' | tab_aggr -k1=1 -k2=7 -s=',' -m=4 > tmp.tab
```
selects entries for the given list of genes in `chromosome.lst` from `cog-24.cog.csv`, sorts entries for each gene in order of the orthology domains within the gene, prints out a (tab-delimited) table where each `gene_id` corresponds to a comma-separated list of `cog_id` to `tmp.tab`.
```
tab_merge chromosome.lst -t= tmp.tab -blank=1
```
prints out the available COG assignments in the chromosome order (adding blanks to column \# 2 for genes without COG assignments).

## COG functional categories in Asgard Archaea
### Source data (https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/):
`cog-24.org.csv` of the following comma-delimited format:
```
genome_id,genome_name,ncbi_taxid,clade
```
`cog-24.cog.csv` of the following comma-delimited format:
```
gene_id,genome_id,...,...,...,...,cog_id,...
```
in columns \## 1, 2 and 7.

`cog-24.def.tab` of the following tab-delimited format:
```
cog_id fun_cat ...
```
in columns \## 1 and 2.
### Commands:
```
tab_select cog-24.org.csv -k1=4 -s1=',' -r='^ASGARD' > tmp.lst
```
prints out a list of ASGARDARCHAEOTA genomes from `cog-24.org.csv` to `tmp.lst`.
```
tab_select cog-24.cog.csv -t= tmp.lst -k1=2 -s1=',' -k2=1 -s2=',' | tab_shuffle -l='1,2,7' -s=',' | sort -u | tab_shuffle -l=3 | tab_merge -t= cog-24.def.tab -k1=1 -k2=1 | tab_regexp -k=2 -r='s/^(.).+/\1/' | tab_aggr -k1=2 | sort
```
selects a subset of COG assignments for these genomes from `cog-24.cog.csv`, trims each line to (tab-delimited) `gene_id genome_id cog_id`, makes this list unique to collapse multiple occurrence of the same orthologous domain within the same gene, extracts `cog_id`, add COG info from `cog-24.def.tab`, trims the functional category (column \# 2) to the single (primary) category, then counts the number of functional category occurrences across all genomes and prints out the sorted list of functional category counts.
