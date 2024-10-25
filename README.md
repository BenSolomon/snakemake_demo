# Installation

- Follow snakemake [installation instructions here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
    - Note snakemake recommends using `mamba` instead of `conda`
- Install the [slurm executor plugin to your snakemake environment](https://pypi.org/project/snakemake-executor-plugin-slurm/)

# Files

- `snakefile_basic.smk`
    - A file containing multiple basic snakemake pipeline operations
    - Uncomment/comment specific sections of the file to demonstrate different operations
    - Includes:
        - Rule syntax anatomy
        - Basic rules
        - Specifying inputs and outputs
        - Chaining rules
        - Using rule all to control outputs
        - Specifying certain output files as temporary
        - Using wildcards
        - Creating a graphic representation of the pipeline
        - Using python scripts
        - Programmatically specifying outputs

- `snakefile_slurm.smk` 
    - A file containing multipe snakemake rules that exemplify the behavior of the SLURM plugin for snakemake
    - Specifically demonstrates the differences between running typical snakemake resource usage **within** an interactive SLURM session vs. using the SLURM plugin to allow snakemake to directly control the SLURM resource management
        - Use within an interactive SLURM session allows for confidently specifying the maximum resources the entire pipeline can use
        - Use of the SLURM executor allows more clear usage of SLURM resources, but it is more challenging to set an overal maximum
    - As above, uncomment/comment specific sections of the file to demonstrate different operations 

- `snakefile_hla.smk`
    - A file containing an example workflow to generate HLA genotypes from fastq files using multiple tools (arcasHLA, optitype, PHLAT)
    - Uses sample fastq files `sample_1.fq.gz` and `sample_2.fq.gz`
    - `hla_bash.sh` represents a SLURM bash script originally used to accomplish a similar workflow (though does have more bells-and-whistles)


# Resources

- [Snakemake docs](https://snakemake.github.io/)
- [Snakemake Slurm plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [UC Berkely Bioinformatics Seminar Snakemake youtube talk](https://www.youtube.com/watch?v=tUTcfoMQl98)
- [YC Lab Snakemake Youtube playlist](https://www.youtube.com/watch?v=Gg0SsEs16Jc&list=PLWhvkMKn3k1zefj7ELcxlukO6AbuP8YCL)