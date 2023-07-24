# VTK2023
# Advanced Technologies for the Life Sciences: Bioinformatics Practical
This markdown document will cover all aspects of the bioinformatics
practical related to the 2022 VTK Advanced Technologies for the Life Sciences.

# Introduction
In this practical we will process 2 RNA sequencing (RNA-seq) datasets
to generate a selection of the analytical results presented in the two
papers associated with the datasets. One of the datasets is a bulk-RNA
dataset while the other is a single cell RNA (scRNA-seq) dataset.

[What is bulk RNA-seq?](https://www.scdiscoveries.com/support/what-is-bulk-rna-sequencing/)

[What is single cell RNA-seq?](https://en.wikipedia.org/wiki/Single-cell_transcriptomics)

## The papers
These are the two papers we will be working with:

- Bulk RNA-seq: [Böstrom et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0188772). Comparative cell cycle
transcriptomics reveals synchronization of developmental
transcription factor networks in cancer cells. PLOS ONE.

- scRNA-seq: [MacParland et al 2018](https://www.nature.com/articles/s41467-018-06318-7). Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nature Communications.

Both papers have been uploaded to the ILLIAS system in `Publications/SequAna`

It would be a good idea to make yourself familiar with the papers before starting the practical work.

## Location
The practical will take place in room M739 from 09.00-16.30 on the 6-8th of December. It is currently unclear whether we will have access to the room on the 9th of December as well.

Ben Hume, the current SequAna bioinformatician will
be running the practical and will be there to assist you.

## Computing setup
You will need access to a computer to complete this practical.
While M739 does contain computers that you could use to login to the
SequAna computational server (where much of the work will occur)
this is a very inconvenient way to work as your home directory is
deleted after every session and it is not possible to install
applications that use Graphical User Interfaces (GUI).

As such, it is strongly recommended that you bring a laptop with
you to complete the practical. There will be 1 spare laptop
available for use during the practical session in M739, but this
laptop cannot leave the room as it belongs to SequAna.

# Objectives

The main objective of this practical is to give you an introduction to the tools that are used by computational biologist/bioinformaticians to generate meaningful results from sequencing data.

The objective of this course is not for you to become proficient or masterful of the techniques we will be covering (we have only 3 or 4 days!), nor to perfectly recreate the figures from the manuscripts. The important part is the journey, not the destination. So take your time. Any proficiency gained in the techniques we cover 
will likely be extremely valuable to you in your career as a research scientist.

To achieve this objective we will work with the sequencing data archived as part of the above mentioned studies to recapitulate several of their key findings.

In doing so we will cover many broad informatic/bioinformatic techniques not limited to:

- Working on the command line interface (CLI)
- Using Conda environments to install programs and packages
- Working with Docker images in Singularity
- Working with core bioinformatic tools to perform:
    - access of archived sequencing data
    - sequencing pre-processing and quality control
    - sequence analysis
- Workflow management with Nextflow
- R scripting to manipulate, analyze and visualize data

I will provide resources for all topics we cover and you are encouraged to look at these
resources if you wish to further your knowledge of the topic.

If you find yourself ahead of the rest of the group, feel free to work on whatever you like or take the time to looking the topics we've covered so far in more detail.

# Structure of the practical
The practical will be divided up by days (1-3). We'll hold the 4th day spare and see how we're getting on. We can be flexible with how or if we use the 4th day.

Each day we will work towards our end goal of recapitulating the results of our chosen studies. But remember, our goal is to learn along the way, not to get to the end. I would rather we take our time on the journey that reach the final figures.

One of the most important skills in computation biology / informatics is the effective
sourcing of reference material. I.e. good googleing!

As such throughout the 3/4 days, while you will be given a structure to follow,
you will also be asked to work out how to do certain tasks on your own.
But don't worry, the SequAna bioinformatician will be there to help you when you get stuck. Much of what you're asked to do will be new to you and may feel challenging - that's totally normal.

If there are any requirements for the day's exercises these will be listed at the beginning of each day's section in a 'requirements' section. For example some packages take a long time to install so it's best to do this setup in advance.

# DAY 1: Installing programs and fetching data

## Requirements
You should already be familiar with working in bash on the command line.
See [here](https://rnabio.org/module-00-setup/0000/08/01/Unix/
) for an excellent resource.

If you don't yet have access to a unix-based operating system you can
use [this](https://cocalc.com/projects?anonymous=terminal).

## Part 1: Connecting to SequAna's computational server: sequana
We will perform much of our analyses on SequAna's
computational server that is creatively named 'sequana'.

The server runs a Linux distribution: Ubuntu.

[What is an operating system?](https://en.wikipedia.org/wiki/Operating_system)

[What is Unix?](https://en.wikipedia.org/wiki/Unix)

[What is Linux?](https://www.linux.com/what-is-linux/)

The vast majority of bioinformatic programs run in Linux,
although many can also be run in Windows or on Mac OS X.

Max OS X is Unix-based and therefore has many similarities to Linux distributions.
The terminal app on Mac OS X offers the user a command line interface (CLI; terminal) very similar to that of Linux distributions.

To connect to sequana we will use ssh.

[What is ssh?](https://www.one.com/en/hosting/what-is-ssh?gclid=Cj0KCQiA1ZGcBhCoARIsAGQ0kkrTFwKizUcQFfXtTDS4WRrrzbZaVqN39hW5ROgMA8ilhZUerT15aWAaAnqlEALw_wcB)

To connect using ssh, we will use ssh-keys.

What are ssh-keys?(https://jumpcloud.com/blog/what-are-ssh-keys)

> **Exercise:** Generate a pair of ssh-keys. Send the public key to Ben for installation on the server along with your username. N.B. if working on Windows, try to use the Open SSH Client on the command line with the command: `ssh-keygen`. For both Windows and mac, the key should be saved to the default location which should be a hidden folder called `.ssh` inside your user's base directory.

> **Exercise:** ssh into the sequana server once Ben has installed your SSH keys. Ask Ben for sequana's IP. (optional) [set up](https://linuxize.com/post/using-the-ssh-config-file/) a 'config' file in '.ssh' so that you don't have to enter your IP and username every time you want to 

## Part 2: Installing software in your user directory
We'll need some programmes if we're going to do some work.

One way to install programs on a Linux system is at the system level.

This allows all users to access the programs. Many of the programs you've already
used are installed at the systems level.

To see where a program is being executed from you can use the command `which`. E.g.:

```
$ which ls
/usr/bin/ls
```

While this might seem convenient, there are many reasons not to install programs at the system level not limited to:
- You must have root access to be able to install the programs
- It is difficult to work with multiple versions of programs (i.e. v0.1.3 vs v0.1.4)

There are several ways to overcome the above issues. One is to download the source code for programs or pre-compiled programs and install these in a local directory. While this solves the problem of root access and the program will only be installed for your use, it still doesn't help us with managing multiple version of programs. It can also be difficult to download the source code of some programs and get it to compile correctly often due to dependency issues. However, sometimes downloading a program and installing it yourself is the only option available to you.

One of the first tasks we will need to undertake for our analysis is getting the sequencing data associated with the studies.

We will download the data from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home).

We will use their [File Downloader Command Line Tool](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html#using-ena-file-downloader-command-line-tool) that is hosted on GitHub [here](https://github.com/enasequence/ena-ftp-downloader/releases).

> **Exercise:** Install the ENA File Downloader in your home directory and make sure you can get it to execute.

We will work with this program later on to download data. But first, let's meet another option for installing programs: [Conda](https://conda.io/projects/conda/en/latest/index.html) - a package manager that allows you to create multiple environments, often specific to a given project or analysis, that can have different packages and programs installed in them.

See [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for managing environmens with conda. Conda can be installed in your home directory without root access (sudo access) and enables you to install programs locally.

> **Exercise:** Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) in your sequana home directory. If you're not sure which options to select, ask! Once it's installed you'll need to start a new ssh session.

Great! We now have the ability to create environments and install programs locally using the 'conda' command.

> **Exercise**: Create a new environment called `kallisto_env` and in it install: 
> - [kallisto](https://anaconda.org/bioconda/kallisto)
>
> Verify that you can run kallisto

We will use this environment later on.
## Part 3: Fetching data, directory permissions and working with symlinks
Now it's time to get some data.

> **Exercise:** Find where the data from each of the studies is deposited. Use the ENA File Downloader that you installed earlier to download data from ONE of the samples from the Böstrom et al 2017 paper.

Wow! Great work, you've now got real sequencing data in your home directory.

I only asked you to get 1 sample's worth of data because it would be redundant to have the same data downloaded multiple times on the server.

And guess what, I've already downloaded the data we need. It's in the following directory on sequana: `/home/VTK/bostrom/raw_reads`.

> **Exercise:** Go take a look at the data. How big is the data?

One of the great things about operating on a Linux system is that its designed to allow concurrent access by multiple users.

> **Exercise**: Take a look at the permissions for the data directory. What does it tell you about access to the directory? What is a group in Linux? Which groups are you a member of?

Now you know where the data is, you can either work directly with that data, or an easier way is to create a shortcut to that data in your own directory. For this purpose, Linux has symbolic links.

[What is a symbolic link (symlink)?](https://www.futurelearn.com/info/courses/linux-for-bioinformatics/0/steps/201767#:~:text=A%20symlink%20is%20a%20symbolic,directory%20in%20any%20file%20system.)

> **Exercise**: Create a directory structure in your home directory to hold the data. Create symlinks to populate the directories with symlinks to the sequencing files.

## Part 4: Some key sequence data formats - fasta, fastq, fastq.gz, sam and bam

There are a few key formats that you should be familiar with in the realms of computational biology.

[What is a fasta file?](https://en.wikipedia.org/wiki/FASTA_format)

[What is a fastq file?](https://en.wikipedia.org/wiki/FASTQ_format)

[What are sam and bam files?](https://www.zymoresearch.com/blogs/blog/what-are-sam-and-bam-files#:~:text=SAM%20files%20are%20a%20type,the%20examples%20for%20this%20section.)

That brings us to the end of the 'first day'. How long did it take us? Hopefully we didn't use up the full day because tomorrow we'll be doing some R in the latter half of the day. So for the remainder of today, I'd like you to spend some time getting familiar with R. This means getting it installed on your system. Either install [Rstudio](https://posit.co/) or [Visual Studio Code](https://code.visualstudio.com/) (my personal favourite). I'm here to help. Once you have that installed R install some of the packages we will be using tomorrow (see the requirement section). Then, if you're not already comfortable with R, or if you're a little rusty, use this time to do a brief R tutorial that covers the basics. R is a fantastic language to get familiar with as a biologist.

# DAY 2: Böstrom et al 2017

## Requirements
Install R on your system or on the sequana server and install BiocManager

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

Install some of the packages we will be using:
- dplyr
- tidyverse
- GenomicFeatures
- RMariaDB
- stringr
- tximport
- DESeq2
- pheatmap
- vsn
- matrixStats

## Part 5: preprocessing:
For Boström there is a single fastq file per sample.

> **Exercise**: Open up the fastq file and have a look at it. Can you deduce the structure of the fastq format from this file? What is the /1 at the end of the reads? How about those quality scores, what use will we make of them?

It's common practice to remove low quality bases at the beginning and ends of reads and to check for and remove adapters?

> **Exercise**: Discuss: Why are there still adapters in my sequences if the sequencing starts after the adapters?

There are a couple of options for removing low quality bases and adapters.

> **Exercise**: Have a look at a couple of the options: [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [fastp](https://github.com/OpenGene/fastp) and [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Specifically, see if you can resolve any differences between the programs. Is one easier to use than the others? Do they all have the same functionality? Let's say you want to remove adapters and trim for quality score at the beginning and end of the read at a phred score of 20.

[What is a phred score?](https://en.wikipedia.org/wiki/Phred_quality_score)

I've worked with each of these over the years. For ease of use and speed, I like fastp.

> **Exercise**: Although we downloaded all of the sequencing files from all of the boström samples, we're only interested in those samples that are HeLa cells. Delete the SymLinks you created that relate to the U20S cells. Do you remeber where to find the information about which SRRs relate to which samples? Do you remember the bostrom_meta.csv file that I created. Perhaps you could look in there also?

> **Exercise**: Install fastp and run it on each of the samples. Set the quality threshold to 20. Think about where the outputs from running fastp will go.

## part 6: Nextflow!
Let's be honest, running fastp once on each of the samples is super annoying. I guess we could generate a bash script to do it for us, but that would be pretty complicated and it would likely run in parallel, and we'd still have to manually take care of the directory structures for the results, and keep track of which versions of the software we used etc. etc. etc. And this was just for one bioinformatic tool! A typical pipeline may have upwards of 10 tools used in sequence, with the outputs combined in certain ways.

What if I told you there was a better way?!

There is! It's called [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html).

> **Exercise**: With Ben's help, write a Nextflow script that will do what we've just done. What do you think? If you were to continue with bioinformatics, would you invest the time to learn Nextflow. What are some of the other benefits of working with Nextflow?

## Part 7: Getting files in and out of the server
Obviously, there will be many occasions where you want to either load files on to the server, or take files off of their server.

We will cover that now.

The most common way to do this is using [scp](https://linuxize.com/post/how-to-use-scp-command-to-securely-transfer-files/) (secure copy).

> **Exercise**: Inspect the output of fastp. Pull down the html file to inspect it.
## Part 8: preparing the reference

> **Exercise**: Discuss: What do you think is the next step after pre-processing? It's important to have a good overview of what you are trying to achieve and how.

It's time to map and count!

When it comes to mapping we have several options available to us.

These can largely be categorized in two:
- map the reads to a genome and them count them
- pseudo-align the reads to a transcriptome

[What is mapping?](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html)

[What is a splice-aware aligner?](https://www.biostars.org/p/175454/)

[What is pseudo-alignment?](https://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html)

> **Exercise**: Discuss: Understanding that pseudo aligners generally map to a transcriptome whereas tradition splice-aware mappers such as STAR map to a genome, can you think of a down side to working with a pseudo-aligner?

We're going to go the pseudo-alignment route with [kallisto](https://pachterlab.github.io/kallisto/): it's fast and does exactly what we need. As we're working with Human, there is a well established set of transcripts to map to (the transcriptome).

I've downloaded the Ensembl transcriptome (release 107) from [here](https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/cdna/). I've downloaded it to here: `/home/VTK/bostrom/reference/`

> **Exercise**: symlink the transcriptome into your directory structure somewhere. Look at the first 10 or so lines of the transciptome using `less` piped into `head`.

We need to do some work to this reference before we use it.

Later on we're going to make a file that maps transcript name to gene name. Don't worry about what that means just now. But what we do need to do is remove the ".1" or ".X" that follows the Ensembl identifier.

We're going to do this using [AWK](https://en.wikipedia.org/wiki/AWK). Even knowing about AWK makes you a next level geek. AWK is sort of a programming language like R or Python. But it acts in a very special way, line by line. The command we're going to run is VERY simple but it gets the job done.

> **Exercise**: In the directory where you have symlinked the reference file to, run the following awk command

```less Homo_sapiens.GRCh38.cdna.all.fa.gz | awk -F '.' '{print $1}' | head```

> What is this doing?
> Now modify the command so that you write the full results to a new file:
>
> ```less Homo_sapiens.GRCh38.cdna.all.fa.gz | awk -F '.' '{print $1}' > Homo_sapiens.GRCh38.cdna.all.short_name.fa```

This new file that we've created `Homo_sapiens.GRCh38.cdna.all.short_name.fa` is the file that we're going to work with for doing our analysis.

> **Exercise:** Gzip the newly created file using `pigz`. If pigz isn't installed, install it.

## Part 9: indexing the reference

Before we map or align, we index our reference.

[Why index?](https://www.biostars.org/p/212594/)

> **Exercise**: Using your kallisto_env environment that you made earlier, index the reference 'short_name' transcriptome.

[What is the GRCh38 assembly?](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19)

## Part 9: pseudo-alignment with kallisto

> **Exercise**: Now perform the pseudo-alignment following the kallisto documentation and incorporating this into your nextflow script. Use 10 cores each time you run kallisto count. You'll need to specify some extra information to kallisto. See if you can figure out what these are. Don't be afraid to have a go running it.

## Part 10: importing your data into R

Once you'ce successfully run kallisto have a look at the output directory of kallisto count. Take a look in some of the files. If you're not familiar with a file's extension, google it.

In abundance.csv you'll see there is a column tpm.

Raw read counts cannot be used to compare expression levels between samples due to the need to account for differences in transcript length, total number of reads per samples, and sequencing biases. Therefore, RNA-seq isoform quantification software summarize transcript expression levels either as TPM (transcript per million), RPKM (reads per kilobase of transcript per million reads mapped), or FPKM (fragments per kilobase of transcript per million reads mapped); all three measures account for sequencing depth and feature length.

> **Exercise**: Research: What is tpm? What are RPKM and FPKM? Do the descriptions make sense?

> **Exercise**: Discuss, what are we going to count? Transcripts or Genes? If we want to count genes, how do we go from transcipts to genes?

There is a package designed specifically to import data from the output of kallisto. It's called [tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html).

We'll use this to import our data. As part of importing the data we will collapse our count data from the transcript level to the gene level. tximport can also help us with that. To do that we need to provide tximport with a dataframe that maps the transcript ID to a gene ID. Several transcripts can originate from a single gene. These are referred to as isoforms. To generate the tx2gene data frame, we can use another package called [makeTxDbFromEnsembl](https://rdrr.io/bioc/GenomicFeatures/man/makeTxDbFromEnsembl.html) from the [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) package. This can all get a little complicated so I'll help you with this.

You will find that there is considerable documentation on the above, both as part of the tximport documentation but also with the [DESEQ2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) program that we will use to do the differential expression analysis.

I've created a file that contains the meta information for the samples here: `/home/VTK/bostrom/bostrom_meta.csv`.

I generated it form [this](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA413699&o=acc_s%3Aa&s=SRR6150369,SRR6150370,SRR6150371,SRR6150372,SRR6150373,SRR6150374,SRR6150375,SRR6150376,SRR6150377,SRR6150378,SRR6150379,SRR6150380,SRR6150381,SRR6150382,SRR6150383,SRR6150384,SRR6150385,SRR6150386,SRR6150387,SRR6150388,SRR6150389). You'll notice that I combined the reads from the same sample.

> **Exercise**: import your outputs from kallisto into R so that we can use DESEQ2 to analyse them. If you're up for the challenge, have a go at doing it yourself following the [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data). Else follow along with me.

Now it's time for us to create the DESeq2 object and perform the analysis.

The DESeq2 documentation is fantastic and I would suggest you open it up now and follow along. To create the first figure from the paper we'll be using some of the standard approaches in this document. OR WILL WE!

> **Exercise**: Have a think about how we create a heat map like the one they're produced in the paper. Have a look at the methods of the paper. Critically appraise them. Do you think they are sufficient to be able to reporduce the findings?
