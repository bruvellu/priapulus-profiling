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

Test run
--------

I will do a test run with Agalma's phylogenetic pipeline using 6 SRA packages.

| species                   | accession | million reads |
| :-----:                   | :-------: | :-----------: |
| _Nanomia bijuga_          | SRX288430 | 52.5          |
| _Apostichopus japonicus_  | SRX122622 | 26.9          |
| _Hormogaster elisae_      | SRX251927 | 30.0          |
| _Mizuhopecten yessoensis_ | SRX220583 | 56.1          |
| _Haemonchus contortus_    | SRX319246 | 29.6          |
| _Penaeus monodon_         | SRX110652 | 29.6          |

```
Metazoa; Eumetazoa; Cnidaria; Hydrozoa; Siphonophora; Physonectae; Agalmatidae; Nanomia; Nanomia bijuga
Metazoa; Eumetazoa; Bilateria; Deuterostomia; Echinodermata; Eleutherozoa; Echinozoa; Holothuroidea; Aspidochirotacea; Aspidochirotida; Stichopodidae; Apostichopus; Apostichopus japonicus
Metazoa; Eumetazoa; Bilateria; Protostomia; Lophotrochozoa; Annelida; Clitellata; Oligochaeta; Haplotaxida; Lumbricina; Hormogastridae; Hormogaster; Hormogaster elisae
Metazoa; Eumetazoa; Bilateria; Protostomia; Lophotrochozoa; Mollusca; Bivalvia; Pteriomorphia; Pectinoida; Pectinoidea; Pectinidae; Mizuhopecten; Mizuhopecten yessoensis
Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Strongylida; Trichostrongyloidea; Haemonchidae; Haemonchinae; Haemonchus; Haemonchus contortus
Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Crustacea; Malacostraca; Eumalacostraca; Eucarida; Decapoda; Dendrobranchiata; Penaeoidea; Penaeidae; Penaeus; Penaeus monodon
```
