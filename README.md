# Advanced Technologies for the Life Sciences (VTK) [BIO-16100] : Bioinformatics Practical
This markdown document will cover all aspects of the bioinformatics
practical related to the VTK Advanced Technologies for the Life Sciences course (BIO-16100).

# Introduction
In this practical, we will process a bulk-RNA sequencing (RNA-seq) dataset
to generate a selection of the analytical results presented in the paper associated with the dataset. single cell RNA (scRNA-seq) is another RNA sequencing approach that will kept as a bonus practice

[What is bulk RNA-seq?](https://www.scdiscoveries.com/support/what-is-bulk-rna-sequencing/)

[What is single-cell RNA-seq?](https://en.wikipedia.org/wiki/Single-cell_transcriptomics)

## The papers
These are the two papers we will be working with, the first one only, while the second one is for your further reading and as a bonus for whoever manages to work on it:

- Bulk RNA-seq: [Böstrom et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0188772). Comparative cell cycle
transcriptomics reveals synchronization of developmental
transcription factor networks in cancer cells. PLOS ONE.

- scRNA-seq: [MacParland et al 2018](https://www.nature.com/articles/s41467-018-06318-7). Single-cell RNA sequencing of the human liver reveals distinct intrahepatic macrophage populations. Nature Communications.

Both papers have been uploaded to the ILLIAS system in `Publications/SequAna`

It would be a good idea to familiarize yourself with the papers before starting the practical work.

## Location
The practical will take place in room M739 from 09.30-16.30 on the 3-6th of December. 

Abdoallah, aka "Abdo," the current SequAna bioinformatician, will run the practical and be there to assist you.

## Computing setup
You'll need access to a computer to complete this practical. While M739 does contain computers still you need to use your laptops. However, we also already set up an online cloud environment (Codespace) for the first part of this practical, but it has limited computing and storage resources. So, we will use SequAna's server for the second part. As such, it is strongly recommended that you bring a laptop with you to complete the practical. 

# Objectives

The main objective of this practical is to introduce you to the tools used by computational biologists/bioinformaticians to generate meaningful results from sequencing data.

The objective of this course is not for you to become proficient or masterful of the techniques we will be covering (we have only 4 or 5 days!), nor to perfectly recreate the figures from the manuscript. The critical part is the journey, not the destination. Please feel free to take your time. Any proficiency gained in the techniques we cover 
will likely be extremely valuable to you as a research scientist.

To achieve this objective, we will use the sequencing data archived as part of the study mentioned above to recapitulate several of their key findings.

In doing so we will cover many broad informatic/bioinformatic techniques not limited to:

- Working on the command line interface (CLI)
- Using Conda environments to install programs and packages
- Working with Docker images in Singularity
- Working with core bioinformatic tools to perform:
    - access to archived sequencing data
    - sequencing pre-processing and quality control
    - sequence analysis
- R scripting to manipulate, analyze, and visualize data
- Workflow management with Nextflow

I will provide resources for all topics we cover and you are encouraged to look at these
resources if you wish to further your knowledge of the topic.

If you find yourself ahead of the rest of the group, you can just work on whatever you like or take the time to look over the topics we've covered so far in more detail.

# Structure of the practical
The practical will be divided into 3-4 days. We'll hold the 5th day spare and see how we're getting on. We can be flexible with how or if we use the 5th day. The 5th day can be used for running the prcticle yourself and preparing the final report, or for the scRNA-seq practical as a bonus session.

Each day, we will work towards our end goal of recapitulating the results of our chosen studies. But remember, our goal is to learn along the way, not to get to the end. I would instead take our time on the journey that reaches the final figures.

One of the most essential skills in computation biology/informatics is the effective
sourcing of reference material. I.e. good googling!

As such throughout the 3/4 days, while you will be given a structure to follow,
you will also be asked to work out how to do certain tasks on your own.
But don't worry, the SequAna bioinformatician will be there to help you when you get stuck. Much of what you're asked to do will be new and may feel challenging - that's normal.

Some packages take a long time to install, so it's best to do this setup in advance. Then, during the first two days, we will install all requirements while we introduce you to the bash and R scripting basics. 

Last but not least, Here are links to tutorial hand-outs for each practical session, Enjoy training.


- DAY 1: [Installing programs bash scripting introduction.](https://github.com/SequAna-Ukon/VTK_CF/wiki/DAY-1:-Installing-programs-bash-scripting-introduction)
- DAY 1: [R scripting (self-studying)](https://github.com/SequAna-Ukon/VTK_CF/wiki/Day-1-(BONUS):-R-scripting)
- DAY 2: [Böstrom et al 2017](https://github.com/SequAna-Ukon/VTK_CF/wiki/DAY-2:-B%C3%B6strom-et-al-2017)
- Day 3: [Functional Enrichment Analysis.](https://github.com/SequAna-Ukon/VTK_CF/wiki/Day-3:-Functional-enrichment-analysis)
- Day 4: [Containers and Workflow management systems.](https://github.com/SequAna-Ukon/VTK_CF/wiki/Day-4:-Containers-and-Workflow-management-systems)
- Day 5 (Bonus): [MacParland et al 2018](https://github.com/SequAna-Ukon/VTK2023/wiki/Day-5:-MacParland-et-al-2018)
