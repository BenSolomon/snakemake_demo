'''
snakemake -s snakefile_basic.smk --use-conda
'''
# # [ ] Rule anatomy
# # This is an example of the types of variables that might
# # be defined in a Snakefile rule. # Don't execute this code.

# rule rule_anatomy:
#     input: 'input file'
#     output: 'output file'
#     conda:
#         'name of conda environment'
#     params:
#         custom_var='custom_value'
#     resources:
#         slots=int,
#         tasks=int
#         mem_mb=int,
#     message: 'Some message to print to STDOUT'
#     log: 'logfile.log' # For a single log
#         stdout='log.stdout',
#         stderr='log.stderr'
#     priority: int # Default 0, any competing rule with a higher priority will be executed first
#     script:
#         'path to python, R, bash, julia, or rust script (maybe others too?)'
#     shell: 
#         '[shell command]'
#     run:
#         [direct python code]

################
# # [ ] Basic rule
# # This will print 'Hello world' to the console

# rule print_hello_world:
#     shell:
#         '''
#         echo 'Hello world'
#         '''

#################
# # [ ] Input and output files
# # This will write 'Hello world' to a file called hello_world.txt

# rule write_hello_world:
#     output:
#         'hello_world.txt'
#     shell:
#         '''
#         echo 'Hello world' > {output}
#         '''

# # Notice that if we attempt to run the Snakefile a second time,
# # `write_hello_world` will not run because the file hello_world.txt already exists.

#################
# # [ ] Chaining rules
# # Example of chaining two rules together

# rule hello_world_A:
#     output:
#         'hello_world_A.txt'
#     shell:
#         '''
#         echo 'Hello world' > {output}
#         '''

# rule hello_world_B:
#     input:
#         'hello_world_A.txt'
#     output:
#         'hello_world_B.txt'
#     shell:
#         '''
#         echo "hello_world_A.txt contains $(cat {input})" > {output}
#         '''

# Note that if this is run as `snakemake -s snakefile_basic.smk` only 
# hello_world_A.txt will be created because rule hello_world_A is assumed to 
# define the output files

# Note that if this is run as `snakemake -s snakefile_basic.smk hello_world_B.txt`,
# snakemake will look at all of the rules to see how hello_world_B.txt can be created
# and recognize that hello_world_A.txt is an input to hello_world_B.txt and run both rules

#################
# # [ ] Rule all
# # Example of using a rule all to specificy the final output files

# rule all:
#     input:
#         'hello_world_B.txt'

# rule hello_world_A:
#     output:
#         'hello_world_A.txt'
#     shell:
#         '''
#         echo 'Hello world' > {output}
#         '''

# rule hello_world_B:
#     input:
#         'hello_world_A.txt'
#     output:
#         'hello_world_B.txt'
#     shell:
#         '''
#         echo "hello_world_A.txt contains $(cat {input})" > {output}
#         '''

# # Now both rules will run if only `snakemake -s snakefile_basic.smk` is run
# # This is because rule all is the first rule and will be used to define all of
# # the files that need to be created. Since it specifies hello_world.B.txt as the
# # file it need in order to run, snakemake will then look at all the other rules
# # to figure out how to make it. 

#################
# # [ ] Temp files
# # Example of removing intermediate files

# rule all:
#     input:
#         'hello_world_final.txt'

# rule hello_world_A:
#     output:
#         temp('hello_world_A.txt')
#     shell:
#         '''
#         echo 'Hello world' > {output}
#         '''

# rule hello_world_B:
#     input:
#         'hello_world_A.txt'
#     output:
#         'hello_world_B.txt'
#     shell:
#         '''
#         echo "Input: {input}" > {output}
#         '''

# rule hello_world_C:
#     input:
#         'hello_world_B.txt'
#     output:
#         temp('hello_world_C.txt')
#     shell:
#         '''
#         echo "Input: {input}" > {output}
#         '''

# rule hello_world_final:
#     input:
#         'hello_world_C.txt'
#     output:
#         'hello_world_final.txt'
#     shell:
#         '''
#         echo "Input: {input}" > {output}
#         '''

# # Note that hello_word_final.txt and hello_world_B.txt are the only file that
# # will be retained because they lack a temp() declaration around their output files

#################
# # [ ] Wildcards
# # We can make this pipeline more flexible by using wildcards
# # This will also allow us to run the same pipeline with multiple inputs

# rule all:
#     input:
#         'testingWildcards1_B.txt',
#         'testingWildcards2_B.txt',

# rule rule_A:
#     output:
#         '{prefix}_A.txt'
#     shell:
#         '''
#         echo 'Hello world' > {output}
#         '''

# rule rule_B:
#     input:
#         '{prefix}_A.txt'
#     output:
#         '{prefix}_B.txt'
#     shell:
#         '''
#         echo "{wildcards.prefix}_A.txt contains $(cat {input})" > {output}
#         '''

#################
# # [ ] Graph representations
# `--rulegraph` will create a directed acyclic graph of how the rules are connected
# `--dag` will create a similar graph but will also show the files that need to be created
# `dot -Tpng > graph.png` can write a visual representation of that graph to file
# dot requires that graphviz is installed

# snakemake -s snakefile_basic.smk --rulegraph | dot -Tpng > rulegraph.png
# snakemake -s snakefile_basic.smk --dag | dot -Tpng > dag.png

#################
# # [ ] Python scripts
# # Directly providing arguments to a python script
# # This rule will run a python script that generates a 
# # file with random numbers and a file with random letters
# # The number of items in each file is defined by 
# # the n_items parameter multiplied by the target 
## file names's suffic number. 
# # E.g numbers_2.txt and n_items=10 should create a file with 20 numbers

# rule all:
#     input:
#         'letters_2.txt',
#         'numbers_2.txt'

# rule python_script:
#     output:
#         'letters_{multiple}.txt',
#         'numbers_{multiple}.txt'
#     params:
#         n_items=10
#     script:
#         'numbers_letters.py'

#################
# # [ ] Expanding output files
# # Expand allows you to specify multiple output files that follow a pattern
# # for example if you have 10 samples numbered sequentially, specify an output 
# # file for all 10 files programmatically using `expand()` 

# # This rule with create a file for each number in range and write
# # that number to the corresponding file

# numbers = range(10)
# rule all:
#     input:
#         expand('{num}.txt', num=numbers)

# rule write_numbers:
#     output:
#         '{num}.txt'
#     run:
#         with open(output[0], 'w') as f:
#             f.write(f'{wildcards.num}\n')

# # Note that this rule also demonstrates how `run:` can be used 
# # to directly run python code in a snakemake rule, rather than 
# # via a script or shell command. 