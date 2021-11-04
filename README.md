# codon_optimization
Script to optimize DNA for recombinant protein expression

To run the codon optimization script (aka reverse_translate.py) from the command line:
python reverse_translate.py FASTAFILE [options: --species (ecoli (default), pichia or saccharomyces) --opt (cai (default) or harmony) --num_outputs (int)]

Examples:
python reverse_translate.py test.fasta
python reverse_translate.py test.fasta --species pichia
python reverse_translate.py test.fasta --opt harmony --num_outputs 100
