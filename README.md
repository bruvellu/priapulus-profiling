bcv-at-brown
============

Records of my visit to the Dunn Lab at Brown University Winter 2014.

Gene expression profiling of priapulid development
--------------------------------------------------

**Objective:** characterize differentially expressed genes during the
development of a priapulid worm using RNAseq data. [Read
more.](priapulid_profiling/description.md)

**Progress:**
- Fixed Fritzen's installation quirks for Agalma.
- Assembled new reference with Agalma 0.3.5.
- Prepared data, calculated read counts, replicate plots.
- Generated STEM profiles.
- Ran differential expression analysis and got significant change in expression
  for many genes that correlate with developmental events.
- Started playing with visualizing the data.

**Plan:**
- Consolidate data analyses and differentially expressed genes.
    - Review RNAseq literature and get to know all the caveats of RNAseq and
      differential expression analyses.
    - Thoroughly understand the statistical strength and weaknesses of the data,
      including all implicit assumptions.
    - Delineate the assumptions of the data and how the differential expression test
      is being executed.
    - Try different DE methods.
- Annotate genes with Gene Ontology and merge with results.
- Select relevant developmental pathways and describe gene activity.
- Explore the possibility of comparing with data sets of other species.
- Think about data visualization.

Phylogenomics of Acoelomorpha
-----------------------------

**Objective:** solve the phylogenetic position of acoels. [Read
more.](acoel_phylogenomics/description.md)

**Progress:**
- Developed a script/module to search and filter SRA records.
- Searched and got a list of SRA packages.
- Tested workflow with a few organisms.
- Defined curation criteria: if unique, yes; if multiple, pick mixed tissues; if
  tissue specific, merge multiple until 40M reads.
- Running assemblies.

**Plan:**
- Continue building assemblies.
- Solve issues like quality offset, no rRNA, too large samples, multiple
  samples.
- Define list of additional RNAseq assemblies.
- Define list of additional genomic assemblies.
- Study phylogenetics.
