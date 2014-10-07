ALDEx2
======
ALDEx tool to examine compositional high-throughput sequence data.

A differential relative count abundance analysis for the comparison of 
two conditions. For example, single-organism and meta-rna-seq 
high-throughput sequencing assays, or of selected and unselected 
values from in-vitro sequence selections. Uses a Dirichlet-multinomial 
model to infer abundance from counts, that has been optimized 
for three or more experimental replicates. Infers sampling 
variation and calculates the expected Benjamini-Hochberg false discovery 
rate given the biological and sampling variation using several parametric
and non-parametric tests. Can to glm and Kruskal-Wallace tests on 
one-way ANOVA style designs.

Note on versioning: ALDEx2_2.0.7.2 was the base for adding the package
to Bioconductor, and the Bioconductor rules indicate that the initial
version must be 0.99.?. Therefore, the current version is 
ALDEx2_0.99.1.tar.gz. This version contains bug fixes and
formatting changes that make it compatible with the Bioconductor 
rules. This version also has parallel execution turned ##OFF# by 
default. To enable parallel execution for aldex.clr, aldex.glm and
aldex.effect set useMC=TRUE when invoking the functions. 

Current version: ALDEx2_0.99.1.tar.gz supplants ALDEx2_2.0.7.2.tar.gz
Current manual: manual/ALDE2_manual.pdf
Current vinette: ALDEx2_vignette.pdf

ALDEx 1.0.4 from Fernandes et al, 2013 PLoS ONE for historical reasons.


