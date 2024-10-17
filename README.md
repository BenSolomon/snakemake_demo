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
        - Using wildcards
        - Creating a graphic representation of the pipeline
        - Using python scripts
        - Programmatically specifying outputs


# Resources

- [Snakemake docs](https://snakemake.github.io/)
- [Snakemake Slurm plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [UC Berkely Bioinformatics Seminar Snakemake youtube talk](https://www.youtube.com/watch?v=tUTcfoMQl98)
- [YC Lab Snakemake Youtube playlist](https://www.youtube.com/watch?v=Gg0SsEs16Jc&list=PLWhvkMKn3k1zefj7ELcxlukO6AbuP8YCL)