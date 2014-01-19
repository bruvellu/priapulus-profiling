bcv-at-brown
============

Records of my visit to the Dunn Lab at Brown University Winter 2014.

Summary
-------

Little is known about the dynamics of gene expression during development in non-model organisms. Next-Generation Sequencing allows a high-throughput sampling of transcripts for any living organism providing a more complete view of gene activity and facilitating comparative studies.

My goal is to analyze the gene expression during the development of a priapulid worm. _Priapulus caudatus_ is a mud-dwelling marine predator known as “penis worms”. Priapulids have a key phylogenetic position belonging to the third main Ecdysozoa lineage together with nematodes and arthropods, the scalidophorans. Abundant fossil record and mitochondrial genes suggest that priapulids have changed little over 500 million years of evolution making them a representative model system for the ancestral ecdysozoan.

Thus, profiling the gene expression of _P. caudatus_ can provide interesting data for a comparative analysis with conventional ecdysozoan model systems such as _Drosophila melanogaster_ and _Caenorhabditis elegans_. 

Questions
---------

* Which patterns of gene expression occur during development? Do they correlate with any developmental event?
* Can these patterns be identified and compared to _D. melanogaster_, for example? What is conserved or missing?

Technical goals
---------------

* Develop hands-on bioinformatics toolset.
* Learn how to process raw data of high-throughput sequencing including stages of Agalma pipeline. What is the rationale and approach behind each step?
* Learn time-series analysis on transcriptome data including statistical tests for identifying differentially expressed genes.

Collection & Sampling
---------------------

We collected adult priapulids by dredging at the Sven Lovén Marine Station in Kristineberg, Sweden. Two large females with ripe gonads were selected for the experiment. We shaked dissected gonads to obtain eggs which we then fertilized using a sperm mix of different males. Cultures were kept at 10 °C.

Every day we picked at least 100 good-looking embryos (usually more) and added directly to a tube with RNAlater.

RNA Extraction
--------------

I extracted the total RNA with phenol/chloroform protocol for female 1 and female 2, separately. Later stages of female 1 were not collected.

|	stage			|	time	|	F1 ng/µL	|	F2 ng/µL	|
|	:----:			|	:---:	|	:-------:	|	:-------:	|
|	oocyte			|	0d		|	153			|	107			|
|	32 cell			|	1d		|	107			|	187			|
|					|	2d		|	306			|	237			|
|	gastrula		|	3d		|	265			|	246			|
|					|	4d		|	124			|	835			|
|	late gastrula	|	5d		|	143			|	215			|
|					|	6d		|	167			|	112			|
|	introvertula	|	7d		|	81			|	303			|
|	pre-hatch larva	|	9d		|	249			|	283			|
|	hatching larva	|	12d		|	-			|	390			|
|					|	15d		|	-			|	485			|
|	lorica larva	|	20d		|	-			|	132			|

Sequencing
----------

### Time-series

Samples of 0, 1, 3, 5, 7, and 9 days for each female were selected for sequencing; 12 samples in total to be sequenced on a single lane. We sent to GeneCore sequencing facility for a Illumina HiSeq2000 Single End with 50 bp reads. Due to low RNA amount, one library preparation failed, female 1 7d. Only 11 samples were sequenced on a single lane.

|	stage			|	time	|	F1 reads	|	F2 reads	|
|	:-----:			|	:----:	|	:--------:	|	:--------:	|
|	oocyte			|	0d		|	40524220	|	19240372	|
|	32 cell			|	1d		|	14403749	|	13153749	|
|	gastrula		|	3d		|	30566620	|	18861835	|
|	late gastrula	|	5d		|	20631875	|	10762813	|
|	introvertula	|	7d		|	-			|	17112976	|
|	pre-hatch larva	|	9d		|	17617125	|	14976016	|

### Reference

RNA samples above were pooled together for a deeper sequencing with a single lane of Illumina HiSeq2000 Paired Ends with 100 bp reads.

|	Paired ends		|	reads								|
|	:---------:		|	:---:								|
|	PcauRNAseq_1	|	[232220195][PcauRNAseq_1_fastqc]	|
|	PcauRNAseq_2	|	[232220195][PcauRNAseq_2_fastqc]	|

[PcauRNAseq_1_fastqc]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Link%20to%20PcauRNAseq_1_fastqc/fastqc_report.html 
[PcauRNAseq_2_fastqc]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Link%20to%20PcauRNAseq_2_fastqc/fastqc_report.html 

Assembly was performed with Agalma 0.3.3 skipping `oases_merge` step.

	agalma transcriptome --id PcauRNAseq --skip oases_merge --outdir /sysdev/s9/bruno/Priapulus_caudatus/RNAseq/agalma/analyses/ > PcauRNAseq.log 2>&1 &