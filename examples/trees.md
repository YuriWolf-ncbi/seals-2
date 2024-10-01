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
