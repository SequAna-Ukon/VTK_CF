# Advanced Technologies for the Life Sciences (VTK2023): Bioinformatics Practical
This markdown document will cover all aspects of the bioinformatics
practical related to the 2023 VTK Advanced Technologies for the Life Sciences.

# Introduction
In this practical, we will process a bulk-RNA sequencing (RNA-seq) dataset
to generate a selection of the analytical results presented in the paper associated with the dataset. single cell RNA (scRNA-seq) is another RNA sequencing approach that will kept as a bonus practice

[What is bulk RNA-seq?](https://www.scdiscoveries.com/support/what-is-bulk-rna-sequencing/)

[What is single-cell RNA-seq?](https://en.wikipedia.org/wiki/Single-cell_transcriptomics)

## The papers
These are the two papers, we will be working with the first one only while the second one is for your further reading and as a bonus for whoever will manage to work on it:

- Bulk RNA-seq: [Böstrom et al 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0188772). Comparative cell cycle
transcriptomics reveals synchronization of developmental
transcription factor networks in cancer cells. PLOS ONE.

- scRNA-seq: [MacParland et al 2018](https://www.nature.com/articles/s41467-018-06318-7). Single-cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations. Nature Communications.

Both papers have been uploaded to the ILLIAS system in `Publications/SequAna`

It would be a good idea to familiarize yourself with the papers before starting the practical work.

## Location
The practical will take place in room M739 from 09.00-16.30 on the 4-8th of December. It is currently unclear whether we will also have access to the room on the 9th of December.

Abdoallah aka "Abdo", the current SequAna bioinformatician will be running the practical and will be there to assist you.

## Computing setup
You'll need access to a computer so you can finish this practical. While M739 does contain computers that use your own laptops. Also, you could use to login to the
SequAna computational server (where much of the work will occur) this is a very inconvenient way to work as your home directory is
deleted after every session and it is not possible to install applications that use Graphical User Interfaces (GUI). However, we also already setted up an online cloud environment (Gitpod) for this practical but it has limited computing and storage resources.

As such, it is strongly recommended that you bring a laptop with
you to complete the practical. 

# Objectives

The main objective of this practical is to introduce you to the tools used by computational biologists/bioinformaticians to generate meaningful results from sequencing data.

The objective of this course is not for you to become proficient or masterful of the techniques we will be covering (we have only 4 or 5 days!), nor to perfectly recreate the figures from the manuscript. The critical part is the journey, not the destination. Please feel free to take your time. Any proficiency gained in the techniques we cover 
will likely be extremely valuable to you as a research scientist.

To achieve this objective we will work with the sequencing data archived as part of the study mentioned above to recapitulate several of their key findings.

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
The practical will be divided up by days (1-5). We'll hold the 5th day spare and see how we're getting on. We can be flexible with how or if we use the 5th day. The 5th day can be used for the scRNA-seq practical as a bonus session.

Each day we will work towards our end goal of recapitulating the results of our chosen studies. But remember, our goal is to learn along the way, not to get to the end. I would instead take our time on the journey that reaches the final figures.

One of the most essential skills in computation biology/informatics is the effective
sourcing of reference material. I.e. good googling!

As such throughout the 3/4 days, while you will be given a structure to follow,
you will also be asked to work out how to do certain tasks on your own.
But don't worry, the SequAna bioinformatician will be there to help you when you get stuck. Much of what you're asked to do will be new to you and may feel challenging - that's totally normal.

Some packages take a long time to install so it's best to do this setup in advance. Then, we will install all requirements during the first two days while we introduce you to the bash and R scripting basics. 

Last but not least, Here are links to tutorial hand-outs for each practical session, Enjoy training.


- DAY 1: [Installing programs bash scripting introduction.](https://github.com/SequAna-Ukon/VTK2023/wiki/DAY-1:-Installing-programs-bash-scripting-introduction)
- DAY 2: [R scripting.](https://github.com/SequAna-Ukon/VTK2023/wiki/Day-2:-R-scripting)
- DAY 3: [Böstrom et al 2017](https://github.com/SequAna-Ukon/VTK2023/wiki/DAY-3:-B%C3%B6strom-et-al-2017)
- Day 4: [Containers and Workflow management systems.](https://github.com/SequAna-Ukon/VTK2023/wiki/Day-4:-Containers-and-Workflow-management-systems)
- Day 5: [MacParland et al 2018(bonus)](https://github.com/SequAna-Ukon/VTK2023/wiki/Day-5:-MacParland-et-al-2018)
