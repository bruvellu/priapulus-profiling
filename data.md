# Collection & Sampling

We collected adult priapulids by dredging at the Sven Lovén Marine Station in
Kristineberg, Sweden. Two large females with ripe gonads were selected for the
experiment. We shook the dissected gonads to obtain eggs which we then
fertilized using a sperm mix of different males. Cultures were kept at 10 °C.
For every day of development, we picked at least 1000 good-looking embryos
(usually more) and added directly to a tube with RNAlater. In total we
collected 9 samples for each female (0-9d, except 8d). We also manage to
collect three additional samples of late larval stages of F2.

Female 1 (F1) was spawned at the station in Sweden. Female 2 (F2) was brought
to the Sars Centre in Bergen, Norway and spawned at the same conditions. 

# RNA Extraction

I extracted the total RNA of each sample using the TriReagent protocol.

| stage           | time  | F1 ng/µL  | F2 ng/µL  |
| :----:          | :---: | :-------: | :-------: |
| oocyte          | 0d    | 153       | 107       |
| 32 cell         | 1d    | 107       | 187       |
|                 | 2d    | 306       | 237       |
| gastrula        | 3d    | 265       | 246       |
|                 | 4d    | 124       | 835       |
| late gastrula   | 5d    | 143       | 215       |
|                 | 6d    | 167       | 112       |
| introvertula    | 7d    | 81        | 303       |
| pre-hatch larva | 9d    | 249       | 283       |
| hatching larva  | 12d   | -         | 390       |
|                 | 15d   | -         | 485       |
| lorica larva    | 20d   | -         | 132       |

# RNA sequencing

The total RNA extractions were shipped to the GeneCore (EMBL Genomics Core
Facilities) for sequencing as detailed below.

## Mixed developmental stages for reference transcriptome

An aliquot of every RNA sample above was pooled together for the deep
sequencing of a reference transcriptome. We used a full lane of Illumina
HiSeq2000 to sequence 100bp paired end reads.

| Paired ends  | reads                 |
| :---------:  | :---:                 |
| Pcau_1       | [232,220,195][Pcau_1] |
| Pcau_2       | [232,220,195][Pcau_2] |

[Pcau_1]: https://rawgit.com/nelas/priapulus-profiling/master/rnaseq/Pcau_1_fastqc.html
[Pcau_2]: https://rawgit.com/nelas/priapulus-profiling/master/rnaseq/Pcau_2_fastqc.html

## Individual developmental stages for expression quantification

For each female, we selected the samples of 0, 1, 3, 5, 7, and 9 days to
quantify the gene expression dynamics. Due to low RNA amount, one library
preparation failed: F1 7d. We thus sequenced 50bp single end reads for the
remaining 11 samples on a single lane of a Illumina HiSeq2000.

| stage           | time | F1 single end reads  | F2 single end reads  |
| :---:           | :--: | :------:             | :------:             |
| oocyte          | 0d   | [40,524,220][Pc1_0d] | [19,240,372][Pc2_0d] |
| 32 cell         | 1d   | [14,403,749][Pc1_1d] | [13,153,749][Pc2_1d] |
| gastrula        | 3d   | [30,566,620][Pc1_3d] | [18,861,835][Pc2_3d] |
| late gastrula   | 5d   | [20,631,875][Pc1_5d] | [10,762,813][Pc2_5d] |
| introvertula    | 7d   | -                    | [17,112,976][Pc2_7d] |
| pre-hatch larva | 9d   | [17,617,125][Pc1_9d] | [14,976,016][Pc2_9d] |

[Pc1_0d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc1_0d_fastqc/fastqc_report.html
[Pc1_1d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc1_1d_fastqc/fastqc_report.html
[Pc1_3d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc1_3d_fastqc/fastqc_report.html
[Pc1_5d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc1_5d_fastqc/fastqc_report.html
[Pc1_9d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc1_9d_fastqc/fastqc_report.html
[Pc2_0d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc2_0d_fastqc/fastqc_report.html
[Pc2_1d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc2_1d_fastqc/fastqc_report.html
[Pc2_3d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc2_3d_fastqc/fastqc_report.html
[Pc2_5d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc2_5d_fastqc/fastqc_report.html
[Pc2_7d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc2_7d_fastqc/fastqc_report.html
[Pc2_9d]: https://dl.dropboxusercontent.com/u/203439/priapulus_caudatus/Pc2_9d_fastqc/fastqc_report.html

