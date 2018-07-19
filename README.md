[![DOI](https://zenodo.org/badge/45686296.svg)](https://zenodo.org/badge/latestdoi/45686296)

HmmUFOtu introduction
=====================
HmmUFOtu is an HMM based Ultra-fast OTU assignment tool for baterial 16S and amplicon sequencing research,
it has two core algorithms, the CSFM-index powered banded-HMM algorithm,
and SEP local phylogenetic-placement based taxonomy assignment algorithm.

The main program 'hmmufotu' takes single or paired-end NGS FASTA/FASTQ reads and generate tab-delimited outputs of the taxonomy assignment results of every read. This program supports native multi-threading and SSE2 and up Vectorization.

HmmUFOtu supports all major DNA substitution models, including
* GTR
* TN93
* HKY85
* F81
* K80
* JC69

HmmUFOtu supports variable mutation rate phylogenetic evaluation using the Discrete Gamma (dΓ) models (Yang 1994).

Citations
---------
> Zheng, Q., Bartow-McKenney, C., Meisel, J. S., & Grice, E. A. (2018). HmmUFOtu: An HMM and phylogenetic placement based ultra-fast taxonomic assignment and OTU picking tool for microbiome amplicon sequencing studies. Genome biology, 19(1), 82. [PubMed ID:     29950165](https://www.ncbi.nlm.nih.gov/pubmed/29950165)

Implementation
--------------
HmmUFOtu is written in pure C++98, and built with the GNU Autotools (autoconfig/automake), and can be easily installed under Linux, Windows and Mac OS X.

Download
--------
You can download the latest release from GitHub at: https://github.com/Grice-Lab/HmmUFOtu/releases.
You can clone or fork and pull the source codes from GitHub at: https://github.com/Grice-Lab/HmmUFOtu.

Dependencies
------------
HmmUFOtu depends on the popular head-only C++ libraries Boost and Eigen3. They are available and often pre-installed on most Linux distributions, and can be easily installed on Windows and Mac OS X.
The ZLIB and Boost-IOSTREAMS libraries are optionally dependent for handling GZIP/GZIP2 compressed files, but are not required and can be disabled.
The JSONCPP C++ library is optionally dependent for formatting HmmUFOtu's assignment output into standard .jplace output; if "jsoncpp" library is found or specified during the configuring step, an optional program **hmmufotu-jplace** will be installed.

Installation
------------
1. Configure installation, by running the command
```bash
./configure
```
You may consider providing additional options, such as `--prefix`, `--exec-prefix`,
`--with-zlib`, `--with-boost`, `--with-libjsoncpp`, etc.

2. Compile and link, by running the command
```bash
make
```
Look for errors and try to resolve them yourself first (by Google).
Contact us only if you are sure it is a bug in our programs.

3. (Optionally) Test, by running the command
```bash
make check
```
It may take a while depending on your processor's speed.

4. Install
```bash
make install
```
You may need root privilege to do it, such as using `sudo`.

Output
------
The main program 'hmmufotu' generates tab-delimited tables (TSV files), and is self explanatory.
One other major program 'hmmufotu-sum' generates TSV format OTU tables (Operational Taxonomic Tables), which is compatitable with 3rd party tools such as QIIME.

Pre-built databases
-------------------
You need to build an HmmUFOtu database before assigning taxonomies to your 16S or other target-loci sequencing reads. You can build your own database using `hmmufotu-build`, or alternatively [download the pre-built databases](https://www.med.upenn.edu/gricelab/hmmufotu.html#databases "Pre-built databases").

Core programs
-------------
Core programs are fundamental tools for taxonomy assignment analysis of 16S and other target-amplicon sequencing data.
Core programs include:
* **hmmufotu-build**		build an HmmUFOtu database with indexed multiple-sequence alignment (MSA), trained HMM profile, and pre-evaluated phylogenetic tree from reference MSA and tree files
* **hmmufotu**			perform HMM-alignment, phylogenetic-placement based taxonomy assignment for single or paired-end NGS reads
* **hmmufotu-sum**		summarize and generate phylogeny-based OTUs and consensus/prior based OTU representatives by summarizing over multiple assignment results (samples)
* **hmmufotu-inspect**		inspect an HmmUFOtu database, and optionally export its contents

Model training programs
-----------------------
HmmUFOtu distributions contain pre-trained DNA substitution models and Dirichlet density/mixture models required for buiding HmmUFOtu databases.
However, the users can train their own models with customized datasets.
Model training programs include:
* **hmmufotu-train-dm**		train an HmmUFOtu prior model using Dirichlet Density/Mixture models with customized data
* **hmmufotu-train-hmm**	train a Banded-HMM model with customized data
* **hmmufotu-train-sm**		train a DNA Substitution Model with customized data

Utility programs
-------------------------
Beside its core functionality, HmmUFOtu can perform many additional analysis using the utility programs.
Utility programs include:
* **hmmufotu-anneal**		anneal primer sequences to an HmmUFOtu database and evaluate the primer efficiency
* **hmmufotu-sim**		generate simulated single or paired-end NGS reads, aligned or un-aligned, using a pre-built HmmUFOtu database
* **hmmufotu-subset**		subset (subsample) an OTUTable so every sample contains the same mimimum required reads, and prune the samples and OTUs if necessary
* **hmmufotu-norm**		normalize an OTUTable so every sample contains the same number of reads, you can generate a relative abundance OTUTable using a constant of 1
* **hmmufotu-merge**  merge two or more OTUTables, redundant OTUs and samples will be aggregated, an optional merged OTU-tree can also be generated providing the corresponding database
* **hmmufotu-jplace**  format HmmUFOtu's assignment output into standard .jplace file for compatibility of third party tools

