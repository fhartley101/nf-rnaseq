# nf-rnaseq

## Introduction

**nf-rnaseq** is a bioinformatics pipeline that can be used to process RNA sequencing data.

The pipeline is built using the [Nextflow](https://www.nextflow.io) Workflow Management System (WfMS), and uses Docker/Singularity containers for enhanced reproducibility. Being written in the [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) language, this pipeline supports the definition of one container per process, thus making it easy to maintain and update software dependencies.

## Pipeline Summary

1. Read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Pseudo-alignment and transcript-level quantification ([Salmon](https://combine-lab.github.io/salmon/))
3. Gene-level summarization ([tximport](https://bioconductor.org/packages/tximport/))
4. Summarize QC Reports ([MultiQC)](https://multiqc.info/)

## Install

1.  Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.6`).
If the [`Conda`](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) 
package manager is already installed in your system, you can use the following 
commands:

    ``` bash
    #Install nextflow
    conda install nextflow

    #Update to the latest version
    nextflow self-update
    ```

2.  Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://docs.sylabs.io/guides/3.11/user-guide/quick_start.html#quick-installation-steps), [`Conda`](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

3. Download the pipeline

    ``` bash
    git clone https://github.com/alebarberis/nf-rnaseq
    ```

4. Start running your own analysis

    ```sh
    nextflow run nf-rnaseq --input <DIR> --outdir <DIR> --salmon_index <DIR> --transcriptome <FILE> -profile <conda/docker/singularity>
    ```

## Usage

``` bash
 nextflow run main.nf --help
```

``` console
=========================================
nf-rnaseq: a nextflow rna-seq pipeline
=========================================
version 0.0.9000

Usage:
The typical command for running the pipeline is as follows:

    nextflow run main.nf --input /path/to/samples --transcriptome /path/to/transcriptome --salmon_index /path/to/index [OPTIONS]

Mandatory arguments:
    --input             DIRPATH             Folder containing FASTQ files, or file with
                                            4 columns: id, read_1, read_2, library_type
    --transcriptome     FILEPATH            FASTA file containing the transcriptome (can be a gzip file)
    --salmon_index      DIRPATH             Folder containing the index on the transcriptome. If empty
                                            a new index will be automatically generated
    --modules           STRING              The pipeline modules to run (default: 'fastqc,quant,multiqc').
                                            Available modules are: fastqc, quant, multiqc

Optional arguments:
    --filext            STRING              Extension of input files (default: fq.gz)
    --suffix1           STRING              Suffix of first file in paired reads (default: _1)
    --suffix2           STRING              Suffix of second file in paired reads (default: _2)
    --concatenate       BOOLEAN             Whether to concatenate input files when multiple files 
                                            per sample id are found (e.g., files from different
                                            lanes)
    --prefix            STRING              Regular expression used to identify groups of multiple
                                            files to concatenate (e.g., --prefix LANE(\d+)_)
    --species           STRING              Species of the samples (e.g., --species hsapiens). 
                                            This parameter is used to create the output sub-folders 
                                            and to download genome/transcriptome data (if required)
    --refdir            DIRPATH             Folder with reference transcriptome and (optional) genome
    --decoys           [FILEPATH]           File containing a set of decoy sequences. If the parameter is
                                            provided without value (i.e., --decoys), a set of decoys
                                            is attempted to be computed from the transcriptome and genome
                                            files
    --genome            FILEPATH            FASTA file containing the genome (can be a gzip file)
    --gtf               FILEPATH            Gene Transfer Format file (can be used to generate a genemap)
    --genemap          [FILEPATH]           File containing a mapping of transcripts to genes. If the 
                                            parameter is provided without a value (i.e., --genemap),
                                            and a GTF file is provided in input, a mapping is attempted 
    --multiqc_config    FILEPATH            Config yaml file for MultiQC
    --outdir            DIRPATH             Output directory (default: ./results)
    --cachedir          DIRPATH             Provide a centralised cache directory for containers (default: ./work)
    --verbose                               Whether to report extra information on progress
    --help                                  Print this usage statement
    --max_cpus          STRING              Maximum amount of allowed cpus (default: 7)
    --max_memory        STRING              Maximum amount of allowed memory (default: '30.GB')
    --max_time          STRING              Maximum amount of execution time (default: '48.h')
```

## Configuration
Nextflow has different configuration sources:

1. Parameters specified on the command line or using a file
2. Configuration files (Nextflow looks for configuration files in multiple locations)
3. Values defined within the pipeline script itself

The sources are ranked in order to decide which settings to apply and avoid conflicts (see [Nextflow's documentation](https://www.nextflow.io/docs/latest/config.html) for further information).

**nf-rnaseq** has the main configuration file (named `nextflow.config`) in the workflow project directory. Other configuration files are located in the `/config` directory: for example, `/config/modules.config` contains DSL2 per module options and publishing paths, while `/config/resources.config` contains resources settings (file modified from the [base.config](https://github.com/nf-core/rnaseq/blob/master/conf/base.config) of the [nf-core rnaseq pipeline](https://nf-co.re/rnaseq).

## Containers
A Container can be seen as a minimal virtual environment or, in simpler words, as a software package containing all the tools needed for a specific task, such as the processing of RNA-seq raw data. The main advantage of containerized software is that it allows the execution of the same analysis on different machines being sure we have the same versions of the computational tools, thus boosting reproducibility. Containers can run on any platform that supports a container runtime.

Being written in the [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) language, this pipeline support the definition of one container per process, thus making it easy to maintain and update software dependencies. The current pipeline default configuration includes profiles (i.e., sets of configuration attributes) for [Conda](https://anaconda.org/), [Docker](https://www.docker.com/), and [Singularity](https://sylabs.io/) containers. A profile can be activated when launching the pipeline execution by using the `-profile` command line option.

``` bash
nextflow run nf-rnaseq --input <DIR> --outdir <DIR> --salmon_index <DIR> --transcriptome <FILE> -profile <conda/docker/singularity>
```

## Custom configuration

### Default parameters and profiles
The default parameters of the pipeline and the defined profiles can be changed by modifying the `/nextflow.config` configuration file.

### Modules
To modify the default configuration of a specific module, you need to identify the module in the `/config/modules.config` configuration file (if missing, create a new entry in the `process` configuration scope) and then simply add/change the container definition. For example, you can change the container for the `FASTQC` module as reported below.

#### Conda
Firstly, check on [Conda](https://anaconda.org/bioconda/fastqc/files) the available versions of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). For example, let's select version 0.11.9. Then, you can modify the configuration accordingly:

``` nextflow
process {
    withName: FASTQC {
        conda = 'bioconda::fastqc=0.11.9'
    }
}
```

#### Docker
Firstly, check on [Quay.io](https://quay.io/repository/biocontainers/fastqc?tab=tags) the available versions of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). For example, let's select version 0.11.9. Then, you can modify the configuration accordingly:

``` nextflow
process {
    withName: FASTQC {
        container = 'quay.io/biocontainers/fastqc:fastqc=0.11.9--hdfd78af_1'
    }
}
```
#### Singularity
Firstly, check on the [Galaxy project](https://depot.galaxyproject.org/singularity) website the available versions of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). The Galaxy project provides all Bioinformatics software from the BioContainers initiative as Singularity prebuilt images. For example, let's select version 0.11.9. Then, you can modify the configuration accordingly:

``` nextflow
process {
    withName: FASTQC {
        container = 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1'
    }
}
```

### Resources
The maximum amount of available memory, CPUs, or allowed execution time, can be modified by changing the related parameters (`max_memory`, `max_cpus`, and `max_time`, respectively) in the `/nextflow.config` configuration file. Process-specific resource requirements can be changed by modifying the `/resources.config` configuration file.

## Credits
The design and implementation of this pipeline was based on the [nf-core rnaseq pipeline](https://nf-co.re/rnaseq) developed by the [nf-core](https://nf-co.re/) project.

The pipeline was written and is maintained by Alessandro Barberis ([@alebarberis](https://github.com/alebarberis)).

## Acknowledgement
Firstly, I would like to express my sincere gratitude to [Prostate Cancer UK](https://prostatecanceruk.org/) for their generous funding, which made it possible for me to develop the first version of this pipeline.
Secondly, I would like to acknowledge and thank Professor [Valentine Macaulay](https://www.nds.ox.ac.uk/team/valentine-macaulay) (head of the [IGF Group](https://www.nds.ox.ac.uk/research/igf-group)) and Professor [Francesca Buffa](https://didattica.unibocconi.it/docenti/cv.php?rif=280949&cognome=BUFFA&nome=FRANCESCA) (head of the [Computational Biology and Integrative Genomics](https://www.oncology.ox.ac.uk/research/research-group/computational-biology-and-integrative-genomics) group) for the invaluable discussions and unwavering support provided throughout the project. And finally, I would like to thank the [Nuffield Department of Surgical Sciences](https://www.nds.ox.ac.uk/) for their consistent assistance.