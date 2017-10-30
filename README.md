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
Some additional features of HmmUFOtu depend on ZLIB and Boost-IOSTREAMS libraries (for handling .gz|.bz2 files), but is not required and can be disabled.

Installation
------------
1. Configure installation, by running the command
```bash
./configure
```
You may consider providing additional options, such as `--prefix`, `--exec-prefix`,
`--with-zlib`, `--with-boost`, etc.

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
You need to build an HmmUFOtu database before assigning taxonomies to your 16S or other target-loci sequencing reads. You can build your own database using `hmmufotu-build`, or alternatively download the pre-built databases below.
* [gg_97_otus_GTR](https://upenn.box.com/shared/static/o146rpg53ebmn3pxikf7zm1uwatez6sl.zip "GreenGenes 97% OUT + GTR")    GreenGenes (v13.8) 97% OTU reference + GTR DNA model
* [gg_97_otus_TN93](https://upenn.box.com/shared/static/ergmvb9gce1t5y4zwrb4soe1kls7it5c.zip "GreenGenes 97% OUT + TN93")    GreenGenes (v13.8) 97% OTU reference + TN93 DNA model
* [gg_97_otus_HKY85](https://upenn.box.com/shared/static/4vdigyt6oywm4lyl0h60f68lc8ag09nt.zip "GreenGenes 97% OUT + HKY85")    GreenGenes (v13.8) 97% OTU reference + HKY85 DNA model

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

