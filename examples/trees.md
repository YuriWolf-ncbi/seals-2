# seals-2 working with trees
## Tidying up a tree and adding species information
### Source data:
`protfam.tre` a Newick format tree file with bare protein IDs as leaf labels.

`prot_species.tab` of the following tab-delimited format:
```
prot_id species
```
(`species` column must be "tree-safe", i.e. not contain spaces, commas and other characters reserved in Newick format).
### Commands:
```
tree_nice protfam.tre > tmp.tre
```
places the bottom trifurcation closest to the geometric center of the tree and sorts tree branches, placing richer clades first.
```
tab_shuffle prot_species.tab -l='1,2,1' -f='%s\t%s|%s' > tmp.tab
```
creates a tab-delimited table containing `prot_id species|prot_id`.
```
tree_rename tmp.tre -t= tmp.tab
```
prints out the tree in Newick format with species names, prefixed to the leaves.
## Parsing a gene tree into low-paralogy subtrees
### Source data:
`genfam.tre` a Newick format tree file with gene IDs as leaf labels.

`gene_genome.tab` of the following tab-delimited format:
```
gene_id genome_id
```
(`genome_id` column must be "tree-safe", i.e. not contain spaces, commas and other characters reserved in Newick format).
### Commands:
```
tab_shuffle gene_genome.tab -l='1,2,1' -f='%s\t%s|%s' > tmp.tab
```
creates a tab-delimited table containing `gene_id genome_id|gene_id`.
```
tree_rename genfam.tre -t= tmp.tab > tmp.tre
```
prints out the tree in Newick format with genome IDs, prefixed to the leaves.
```
tree_dismember tmp.tre -w=1 -s='\|' -n= subfam > subfam.lst
```
parses `tmp.tre` into high-coverage low-paralagy subtrees, named `[[subfam.1.tre], subfam.2.tre, ...], [subfam.xxx.tre]` and prints out the list of subtrees into `subfam.lst`.
