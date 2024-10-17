# Script creates one list of n-random letters and 
# another list of n-random numbers then writes both to 
# separate files. Can be run from command line or snakemake.

import argparse
import random
import string

# Argument parsing for command line
def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_items', type=int, default=100)
    parser.add_argument('--letters_path', type=str, default='letters.txt')
    parser.add_argument('--numbers_path', type=str, default='numbers.txt')
    return parser.parse_args()
    # python numbers_letters.py --n_items 1000 --letters_path letters.txt --numbers_path numbers.txt

def main():
    # Argument parsing from snakemake
    if 'snakemake' in globals():
        # Parsing inputs and outputs from snakemake
        letters_path= snakemake.output[0] # type: ignore
        numbers_path= snakemake.output[1] # type: ignore
        # Parsing parameters from snakemake
        n_items = snakemake.params.n_items  # type: ignore
        # Parsing wildcards from snakemake
        multiple = int(snakemake.wildcards[0]) # type: ignore
        n_items = n_items * multiple
    
    # Argument parsing from command line    
    else:
        args = arg_parser()
        n_items = args.n_items
        letters_path = args.letters_path
        numbers_path = args.numbers_path
    
    # Generate random letters and numbers    
    letters = [random.choice(string.ascii_letters) for _ in range(n_items)]
    numbers = [random.randint(0, 9) for _ in range(n_items)]
    
    # Write numbers and letters to files   
    with open(letters_path, 'w') as f:
        for letter in letters:
            f.write(letter + '\n')    
    with open(numbers_path, 'w') as f:
        for number in numbers:
            f.write(str(number) + '\n')

if __name__ == '__main__':
    main()