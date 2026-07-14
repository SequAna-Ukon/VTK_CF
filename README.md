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

The main objective of this practical course is to introduce you to the tools and computational approaches used by bioinformaticians to extract meaningful biological insights from high-throughput sequencing data.

The aim of this course is not for you to become an expert in the techniques we cover—we only have four to five days—nor is it to perfectly reproduce every figure from the selected publications. Instead, the focus is on understanding the analysis process and developing the practical skills needed to perform each step independently. The journey is more important than the destination, so take your time and focus on learning the concepts and reasoning behind each analysis. Any proficiency you gain in these techniques will be highly valuable throughout your research career.

To achieve these objectives, we will use publicly available sequencing data from the selected studies to reproduce several of their key findings while learning the underlying computational methods.

Throughout the practical, we will cover a range of core computational biology and bioinformatics topics, including:

- Working in the Linux command-line interface (CLI)
- Using Conda environments to install and manage software packages
- Working with Docker containers through Singularity/Apptainer
- Using core bioinformatics tools for:
    - Accessing publicly archived sequencing data
    - Sequence data preprocessing and quality control
    - RNA-seq data analysis
- R scripting for data manipulation, statistical analysis, and visualization
- Workflow management using Nextflow

Resources and reference materials will be provided for every topic covered during the course. You are encouraged to explore these materials further if you wish to deepen your understanding beyond the practical sessions.

If you finish a session ahead of the rest of the group, feel free to continue exploring the provided resources, revisit previous exercises in greater depth, or work on your own analyses or research questions.

# Structure of the practical
The practical course will be divided into 4 core training days, with Day 5 reserved as a flexible project day. This final day will allow everyone to complete their own RNA-seq analysis project on cell apoptosis, finalize any remaining analyses, prepare figures, and work on their final report. We can be flexible with how this day is used depending on the progress of the class.

Each day, we will work towards our overall goal of reproducing the key findings from selected published studies. However, the primary objective of this course is not simply to reach the final results, but to understand the analysis process and develop the computational skills required to perform each step independently. We will therefore take our time throughout the journey, focusing on the methods and reasoning behind each analysis.

One of the most important skills in computational biology and bioinformatics is the ability to efficiently find and evaluate reference material—that is, developing effective literature searching and "Googling" skills. Throughout the course, you will be provided with a structured guide, but you will also be encouraged to discover how to perform certain tasks on your own. Don't worry if you get stuck—the SequAna bioinformatics team will be available throughout the practical sessions to provide guidance and support. Much of the material will be new and may seem challenging at first, and that is a normal part of the learning process.

Some software packages require considerable time to install and configure, so we recommend completing as much of the setup as possible before the course begins. During the first two days, we will ensure that everyone has the required software installed while introducing the fundamentals of Bash and R scripting.

Finally, below are the tutorial handouts for each practical session. We hope you enjoy the training!


- **DAY 1:** [Installing programs bash scripting introduction.](https://github.com/SequAna-Ukon/VTK_CF/wiki/DAY-1:-Installing-programs-bash-scripting-introduction)
- **DAY 1:** [R scripting (self-studying)](https://github.com/SequAna-Ukon/VTK_CF/wiki/Day-1-(BONUS):-R-scripting)
- **DAY 2:** [Böstrom et al 2017](https://github.com/SequAna-Ukon/VTK_CF/wiki/DAY-2:-B%C3%B6strom-et-al-2017)
- **Day 3:** [Functional Enrichment Analysis.](https://github.com/SequAna-Ukon/VTK_CF/wiki/Day-3:-Functional-enrichment-analysis)
- **Day 4:** [Containers and Workflow management systems.](https://github.com/SequAna-Ukon/VTK_CF/wiki/Day-4:-Containers-and-Workflow-management-systems)
- **Day 5:** [Student project](https://github.com/SequAna-Ukon/VTK_CF/wiki/Day-5:-Student-project).
- **Bonus (Self-study):** [MacParland et al 2018](https://github.com/SequAna-Ukon/VTK_CF/wiki/Bonus-(Self%E2%80%90study):-MacParland-et-al-2018)
