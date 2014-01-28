RSEM expression analyses
========================

While in Bergen I ran the basic expression analysis for female 2 following [RSEM steps](http://deweylab.biostat.wisc.edu/rsem/).

Prepare reference sequences
---------------------------

RSEM needs to prepare and build bowtie indices for the reference sequences (see [docs](http://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html)). I'm doing this using Agalma's cleaned assembly.

	cd /sysdev/s9/bruno/Priapulus_caudatus/RNAseq_profiling/
	rsem-prepare-reference /sysdev/s9/bruno/Priapulus_caudatus/RNAseq/agalma/analyses/PcauRNAseq/12/assembly_143959988_trinity.clean.fa references/agalma_cleaned > references/agalma_cleaned.log 2>&1 &

Calculate expression
--------------------

Only processing for female 2 samples. Command below was ran for each time-point.

	/usr/local/bin/rsem/rsem-calculate-expression --fragment-length-mean 52.0 --fragment-length-sd 0.0 -p 12 ~/Priapulus_caudatus/RNAseq_profiling/agalma/data/raw/Pc2_0d.txt references/agalma_cleaned expression/Pc2_0d > expression/Pc2_0d.log 2>&1 &

Visualization
-------------

	cd expression

### Generating a wiggle file
	
	/usr/local/bin/rsem/rsem-bam2wig Pc2_0d.transcript.sorted.bam ../wiggle/Pc2_0d.wig Pc2_0d

### Exporting wiggle plot
	
	/usr/local/bin/rsem/rsem-plot-transcript-wiggles Pc2_0d ../wiggle/Pc2_0d.ids.txt ../wiggle/Pc2_0d.pdf

### Make model

	/usr/local/bin/rsem/rsem-plot-model Pc2_0d ../model/Pc2_0d.pdf

Differential expression analysis
--------------------------------

	https://github.com/bli25wisc/RSEM#-differential-expression-analysis

	mkdir dea
	/usr/local/bin/rsem/rsem-generate-ngvector ../references/agalma_cleaned.transcripts.fa agalma_cleaned

	/usr/local/bin/rsem/rsem-generate-data-matrix ../expression/Pc2_0d.genes.results ../expression/Pc2_1d.genes.results ../expression/Pc2_3d.genes.results ../expression/Pc2_5d.genes.results ../expression/Pc2_7d.genes.results ../expression/Pc2_9d.genes.results > Pcau.counts.matrix

	/usr/local/bin/rsem/rsem-run-ebseq --ngvector agalma_cleaned.ngvec Pcau.counts.matrix 1,1,1,1,1,1 Pcau_genes.results

	/usr/local/bin/rsem/rsem-control-fdr Pcau_genes.results 0.05 Pcau_genes.de.txt

Done with analysis, patterns file has the possible patterns, .de.txt has the significant selection below 0.05.
