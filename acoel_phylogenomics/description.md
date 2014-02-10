Phylogenomics of Acoelomorpha
=============================

1. Fetch high quality RNAseq data from SRA.
2. For each unique organism curate packages to be used.
3. Include well annotated genomes.
4. Insert raw data records into Agalma.
5. Download data.
6. Run Agalma phylogenomics pipeline.

Fetching SRA packages
---------------------

I develop the python module [fetch_sra.py](fetch_sra.py) to search and fetch
packages from [SRA Database](http://www.ncbi.nlm.nih.gov/sra/). Using a defined
set of criteria I (1) searched SRA packages, (2) filtered packages that matched
these parameters. The commands are described in [acoel_sra.py](acoel_sra).

Search criteria:

- Library Strategy: RNA-seq
- Platform: Illumina or 454
- Organism: Metazoa, but not Vertebrata or Insects
- Modification data: After year 2000

Post-search filter criteria:

- Library Layout: paired ends (with nreads=2)
- Read Length: greater or equal than 70 bp

My initial search yield 3117 packages and 736 met our post-search filter
requirements. The packages were written to
[sra_paired_gte70bp.csv](sra_paired_gte70bp.csv). These packages cover 91
unique organisms, shown in
[unique_sra_paired_gte70bp.csv](unique_sra_paired_gte70bp.csv).

Package curation
----------------

Due to the large variety of sampled tissues, some manual curation of packages
is needed to selected the most representative sequencing experiments.

TODO
----

- Deal with previous Illumina ascii offset of 64.
- Merge representative samples into single FASTQ pair.
- Import previous Agalma assemblies from other databases.
- Import other assemblies.
- Import Genomes.

Test run
--------

I will do a test run with Agalma's phylogenetic pipeline using 6 SRA packages.

| species                   | accession | run       | million reads | offset |
| :-----:                   | :-------: | :-:       | :-----------: | :----: |
| _Penaeus monodon_         | SRX110652 | SRR388222 | 29.6          | 33     |
| _Onchocerca volvulus_     | ERX200391 | ERR225731 | 39.2          |        |
| _Crassostrea virginica_   | SRX118365 | SRR404226 | 26.4          |        |
| _Hormogaster elisae_      | SRX251927 | SRR786599 | 30.0          | 33     |
| _Apostichopus japonicus_  | SRX122622 | SRR414930 | 26.9          | 33     |
| _Nanomia bijuga_          | SRX288430 | SRR871527 | 52.5          | 33     |
|                           |           |           |               |        |
| _Haemonchus contortus_    | SRX319246 | SRR928056 | 29.6          | 64     |
| _Mizuhopecten yessoensis_ | SRX220583 | SRR653778 | 56.1          | 64     |

```
Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Crustacea; Malacostraca; Eumalacostraca; Eucarida; Decapoda; Dendrobranchiata; Penaeoidea; Penaeidae; Penaeus; Penaeus monodon
Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Nematoda; Chromadorea; Spirurida; Filarioidea; Onchocercidae; Onchocerca; Onchocerca volvulus
Metazoa; Eumetazoa; Bilateria; Protostomia; Lophotrochozoa; Mollusca; Bivalvia; Pteriomorphia; Ostreoida; Ostreoidea; Ostreidae; Crassostrea; Crassostrea virginica
Metazoa; Eumetazoa; Bilateria; Protostomia; Lophotrochozoa; Annelida; Clitellata; Oligochaeta; Haplotaxida; Lumbricina; Hormogastridae; Hormogaster; Hormogaster elisae
Metazoa; Eumetazoa; Bilateria; Deuterostomia; Echinodermata; Eleutherozoa; Echinozoa; Holothuroidea; Aspidochirotacea; Aspidochirotida; Stichopodidae; Apostichopus; Apostichopus japonicus
Metazoa; Eumetazoa; Cnidaria; Hydrozoa; Siphonophora; Physonectae; Agalmatidae; Nanomia; Nanomia bijuga

Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Strongylida; Trichostrongyloidea; Haemonchidae; Haemonchinae; Haemonchus; Haemonchus contortus
Metazoa; Eumetazoa; Bilateria; Protostomia; Lophotrochozoa; Mollusca; Bivalvia; Pteriomorphia; Pectinoida; Pectinoidea; Pectinidae; Mizuhopecten; Mizuhopecten yessoensis
```

```sh
# Set resources.
export BIOLITE_RESOURCES="database=/sysdev/s9/bruno/acoel_phylo/biolite.sqlite,threads=24,memory=800G,outdir=/sysdev/s9/bruno/acoel_phylo/analyses"

cd /sysdev/s9/bruno/acoel_phylo/data/

sra import --clean --email organelas@gmail.com --gzip --metadata SRX110652
sra import --clean --email organelas@gmail.com --gzip --metadata ERX200391
sra import --clean --email organelas@gmail.com --gzip --metadata SRX118365
sra import --clean --email organelas@gmail.com --gzip --metadata SRX251927
sra import --clean --email organelas@gmail.com --gzip --metadata SRX122622
sra import --clean --email organelas@gmail.com --gzip --metadata SRX288430

catalog insert --id SRX110652 --paths SRR388222_1.fastq.gz SRR388222_2.fastq.gz
catalog insert --id ERX200391 --paths ERR225731_1.fastq.gz ERR225731_2.fastq.gz
catalog insert --id SRX118365 --paths SRR404226_1.fastq.gz SRR404226_2.fastq.gz
catalog insert --id SRX251927 --paths SRR786599_1.fastq.gz SRR786599_2.fastq.gz
catalog insert --id SRX122622 --paths SRR414930_1.fastq.gz SRR414930_2.fastq.gz
catalog insert --id SRX288430 --paths SRR871527_1.fastq.gz SRR871527_2.fastq.gz

cd /sysdev/s9/bruno/acoel_phylo/scratch/

agalma transcriptome --id SRX110652 > Pmon.out 2>&1 &
agalma transcriptome --id ERX200391 > Ovol.out 2>&1 &
agalma transcriptome --id SRX118365 > Cvir.out 2>&1 &
agalma transcriptome --id SRX251927 > Heli.out 2>&1 &
agalma transcriptome --id SRX122622 > Ajap.out 2>&1 &
agalma transcriptome --id SRX288430 > Nbij.out 2>&1 &

# Load.
agalma load --id SRX110652 --previous SRX110652
agalma load --id ERX200391 --previous ERX200391
agalma load --id SRX118365 --previous SRX118365
agalma load --id SRX251927 --previous SRX251927
agalma load --id SRX122622 --previous SRX122622
agalma load --id SRX288430 --previous SRX288430

# Homologize
agalma homologize 7 8 9 10 11 12 --restart homologize.chk --id PhyloTest > PhyloTest.out 2>&1 &

agalma multalign --id PhyloTest
agalma genetree --id PhyloTest
agalma treeprune --id PhyloTest
agalma multalign --id PhyloTest
agalma supermatrix --id PhyloTest
agalma genetree --id PhyloTest --raxml_flags="-o Nanomia_bijuga"

# Reports
agalma report --id PhyloTest --outdir reports/PhyloTest
agalma resource_report --id PhyloTest --outdir reports/PhyloTest
agalma phylogeny_report --id PhyloTest --outdir reports/PhyloTest
```
