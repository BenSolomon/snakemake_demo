#################
# # [ ] Processor usage with interactive SLURM
# # When managing resources with Snakemake, command line --cores controls the maximum cores for the pipeline
# # When managing resources with Snakemake, snakemake `threads:` controls the number of cores for a rule

# Start interactive SLURM with :
# # srun --pty --cpus-per-task=4 --mem=20G --time=01:00:00 --job-name=SLURM_INT --partition=khatrilab --nodes=1 bash

# Run snakemake without direct call to SLURM
# # snakemake -s snakefile_slurm.smk --cores 2
# # snakemake -s snakefile_slurm.smk --cores 2 --resources mem_mb=5000

# Monitor with htop

# n=2 # Number of iterations (in this case, jobs)
# t=4 # Number of threads requested by rule, try 1 and 4

# rule all:
#     input:
#         expand(f"slurmIntProc_t{t}_n{n}_{{jobID}}.txt", jobID=range(n))

# rule cores:
#     threads: t
#     output:
#         f"slurmIntProc_t{t}_n{n}_{{jobID}}.txt"
#     shell:
#         '''
#         echo nproc: $(nproc) > {output}
#         '''

# E.g. n=2, t=4 - Snakemake will attempt to run 2 jobs, each with 4 cores. 
# # However, --cores limits max job cores to 2, so only one job will run at a time and it will be limited to the use of 2 cores
# E.g. n=2, t=1 - Snakemake will attempt to run 2 jobs, each with 1 core. 
# # Since --cores max is 4, both jobs will run simultaneously with 1 core (2 total cores in use)

#################
# # [ ] Memory usage with interactive SLURM
# # Within the interactive SLURM session, jobs run with snakemake will not exceede the limits of the SLURM session

# Start interactive SLURM with :
# # srun --pty --cpus-per-task=4 --mem=20G --time=01:00:00 --job-name=SLURM_INT --partition=khatrilab --nodes=1 bash

# Run snakemake without direct call to SLURM
# # snakemake -s snakefile_slurm.smk --cores 2
# # snakemake -s snakefile_slurm.smk --cores 2 --resources mem_mb=5000

# Monitor with htop

# Will use stress-ng to test memory usage (installed to snakemake conda env)
# Need to ensure stress-ng --timeout is long enough to allow for desired allocation size

# n=2 # Number of iterations (in this case, jobs)
# t=1 # Number of threads requested by rule
# m=40 # Test 1, 10, 20, 40

# rule all:
#     input:
#         expand(f"slurmIntMem_t{t}_n{n}_m{m}_{{jobID}}.txt", jobID=range(n))

# rule cores:
#     threads: t
#     output:
#         f"slurmIntMem_t{t}_n{n}_m{m}_{{jobID}}.txt"
#     resources:
#         mem_mb=10000
#     params:
#         mem=f'{m}G'
#     shell:
#         '''
#         echo nproc: $(nproc) > {output}
#         /usr/bin/time -v stress-ng --vm 1 --vm-bytes {params.mem} --timeout 10s 2>&1 | grep "Maximum resident" | awk '{{print "Maximum resident set size: " $6 / 1024 " MB"}}' >> {output}
#         '''

# E.g. Running a 20G memory job (m="20G") when slurm --mem=20G, will run fine, Maximum resident set size: 20480.5 MB
# E.g. Running a 40G memory job (m="40G") when slurm --mem=20G, will fail, Maximum resident set size: 62.4102 MB

# --resources mem_mb and resources:mem_mb do not actually seem to control memory usage, only how jobs are scheduled 
# # A job using --vm-bytes 20G will run fine on a slurm session with --mem=20G even if --resources mem_mb=10.
# # However, with --resources mem_mb=10, a rule with mem_mb=20 will only run one-at-a-time, while a rule with mem_mb=1 could run 10-at-a-time
# # This can create a situation where two tasks are competing for the same maximum memory set by the SLURM session
# However, at least you know that you won't have runawawy memory use that will use up all available SLURM resources

    
#################
# # [ ] SLURM directly from snakemake
# # When managing resources with Snakemake-SLURM, command line --jobs is the rate limiting parameter
# # When managing resources with Snakemake-SLURM, command line --cores controls the max cores for EACH JOB
# # So the total core usage with be --jobs * --cores 

# `threads:` has no effect on resources

# snakemake \
#     -s snakefile_slurm.smk \
#     --executor slurm \
#     --cores 2 \
#     --jobs 10 \
#     --default-resources \
#     slurm_account=khatrilab \
#     slurm_partition=khatrilab 

# Note that --jobs is required for the SLURM executor

# Monitor with:
# watch -n 2 squeue -o %.8i,%.9P,%.8j,%.8u,%.2t,%.4C,%.6m,%.14N,%.10M,%.20S,%.40k -u solomonb

# n=20 # Number of iterations (in this case, jobs)
# t=20 # Number of threads requested by rule

# rule all:
#     input:
#         expand(f"slurmBatchProc_t{t}_n{n}_{{jobID}}.txt", jobID=range(n))

# rule cores:
#     threads: t
#     output:
#         f"slurmBatchProc_t{t}_n{n}_{{jobID}}.txt"
#     shell:
#         '''
#         sleep 5
#         echo "nproc $(nproc)" > {output}
#         echo "SLURM_CPUS_PER_TASK $SLURM_CPUS_PER_TASK" >> {output}
#         echo "SLURM_JOB_CPUS_PER_NODE $SLURM_JOB_CPUS_PER_NODE" >> {output}
#         echo "SLURM_TASKS_PER_NODE $SLURM_TASKS_PER_NODE" >> {output}
#         echo "SLURM_CPUS_ON_NODE $SLURM_CPUS_ON_NODE" >> {output}
#         echo "SLURM_JOB_CPUS_PER_NODE $SLURM_JOB_CPUS_PER_NODE" >> {output}
#         echo "SLURM_NTASKS $SLURM_NTASKS" >> {output}
#         echo "SLURM_NPROCS $SLURM_NPROCS" >> {output}
#         '''

# The executor will choose how many jobs to run based on the --jobs paramater
# Each job will be allowed a max number of cores set by --cores

# E.g. n=20, t=10, --jobs 10 --cores 2 - 
# # Snakemake will attempt to run n=20 jobs (n), but will only run 10 at a time because of --jobs 10
# # Snakemake will assign 2 cores to each job because of --cores 2
# # threads:10 has no impact on resource management

#################
# # [ ] SLURM specific resources flags in snakemake
# # resources:cpus_per_task overrides --cores

# snakemake \
#     -s snakefile_slurm.smk \
#     --executor slurm \
#     --cores 50 \
#     --jobs 4 \
#     --default-resources \
#     slurm_account=khatrilab \
#     slurm_partition=khatrilab-gpu 


# n=4 # Number of iterations (in this case, jobs)
# t=50 # Number of threads requested by rule
# cpt=30 # CPUs per task

# rule all:
#     input:
#         expand(f"slurmBatchProc_t{t}_cpt{cpt}_n{n}_{{jobID}}.txt", jobID=range(n))

# rule cores:
#     threads: t
#     resources:
#         cpus_per_task=cpt,
#         slurm_extra="--nodelist=bmir-ct-gpu-2"
#     output:
#         f"slurmBatchProc_t{t}_cpt{cpt}_n{n}_{{jobID}}.txt"
#     shell:
#         '''
#         sleep 5
#         echo "nproc $(nproc)" > {output}
#         echo "SLURM_CPUS_PER_TASK $SLURM_CPUS_PER_TASK" >> {output}
#         echo "SLURM_JOB_CPUS_PER_NODE $SLURM_JOB_CPUS_PER_NODE" >> {output}
#         echo "SLURM_TASKS_PER_NODE $SLURM_TASKS_PER_NODE" >> {output}
#         echo "SLURM_CPUS_ON_NODE $SLURM_CPUS_ON_NODE" >> {output}
#         echo "SLURM_JOB_CPUS_PER_NODE $SLURM_JOB_CPUS_PER_NODE" >> {output}
#         echo "SLURM_NTASKS $SLURM_NTASKS" >> {output}
#         echo "SLURM_NPROCS $SLURM_NPROCS" >> {output}
#         '''

# Limiting jobs specifically to node bmir-ct-gpu-2 which has 40 cpus
# `threads:` has no effect on resources, since setting t=50 will not cause any resource failure
# `resources:cpus_per_task` controls the number of cores per job, since cpt=50 will cause resource failure
# nproc will nearly match resources:cpus_per_task (may be off by one)

# If  `resources:cpus_per_task` is not specified --cores will control cores per job. It will not overrule `resources:cpus_per_task`
# # resources:cpus_per_task=50 with --cores 1 will cause failure because cpus_per_task exceeds bmir-ct-gpu-2 cpus even though --cores is set to 1
# # --cores 50 without resources:cpus_per_task specified will cause failure because --cores exceeds bmir-ct-gpu-2 cpus
# # resources:cpus_per_task=1 with --cores 50 will not cause failure because cpus_per_task does not exceed bmir-ct-gpu-2 cpus even though --cores is set to 50

#################
# # [ ] SLURM executor memory management

# # Total memory will be either 

# snakemake \
#     -s snakefile_slurm.smk \
#     --executor slurm \
#     --cores 2 \
#     --jobs 4 \
#     --default-resources \
#     slurm_account=khatrilab \
#     slurm_partition=khatrilab-gpu 


# n=4 # Number of iterations (in this case, jobs)
# cpt=4 # CPUs per task
# mpc=1000 # Memory (MB) per CPU


# rule all:
#     input:
#         expand(f"slurmBatchMem_mpc{mpc}_cpt{cpt}_n{n}_{{jobID}}.txt", jobID=range(n))

# rule cores:
#     resources:
#         cpus_per_task=cpt,
#         mem_mb_per_cpu=mpc,
#         slurm_extra="--nodelist=bmir-ct-gpu-3"
#     output:
#         f"slurmBatchMem_mpc{mpc}_cpt{cpt}_n{n}_{{jobID}}.txt"
#     params:
#         mem='3G'# Try 3G and 5G
#     shell:
#         '''
#         echo "nproc $(nproc)" > {output}
#         echo "SLURM_CPUS_PER_TASK $SLURM_CPUS_PER_TASK" >> {output}
#         echo "SLURM_JOB_CPUS_PER_NODE $SLURM_JOB_CPUS_PER_NODE" >> {output}
#         echo "SLURM_TASKS_PER_NODE $SLURM_TASKS_PER_NODE" >> {output}
#         echo "SLURM_CPUS_ON_NODE $SLURM_CPUS_ON_NODE" >> {output}
#         echo "SLURM_JOB_CPUS_PER_NODE $SLURM_JOB_CPUS_PER_NODE" >> {output}
#         echo "SLURM_NTASKS $SLURM_NTASKS" >> {output}
#         echo "SLURM_NPROCS $SLURM_NPROCS" >> {output}
#         echo "SLURM_MEM_PER_CPU $SLURM_MEM_PER_CPU" >> {output}
#         /usr/bin/time -v stress-ng --vm 1 --vm-bytes {params.mem} --timeout 5s 2>&1 | grep "Maximum resident" | awk '{{print "Maximum resident set size: " $6 / 1024 " MB"}}' >> {output}
#         '''

# Unlike with interactive session, attempting to run a process with requiring more memory than 
# available based on the resources set by Snakemake-SLURM will cause the job to fail, rather than
# just scale down the memory usage.

# E.g. If running 4 jobs, each with 4 cores, and requesting 1G RAM per core, should result in 4 jobs with 4 CPUs and 4G RAM each
# If you run a process in those jobs requiring more than 4G RAM, the job will fail 
