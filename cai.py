#!/usr/bin/env python
import random
import sys
import math
import argparse

parser = argparse.ArgumentParser(description='Manage script input')
parser.add_argument("fasta_files", default=None, nargs="*", help="Path to the input xlsx file(s).")
parser.add_argument("-s", "--species", default='ecoli', help="Expression species to optimize for. Currently only accepts ecoli and pichia")
# parser.add_argument("-o", "--outpath", default='.', help="Path to output folder")
args = parser.parse_args()

if len(args.fasta_files) >= 1:
    fasta_files = args.fasta_files
else:
    # files = glob.glob('.fasta')
    sys.exit('You have to provide a fasta file to optimize!\nUsage: python reverse_translate.py FASTA_FILE [optional: --species pichia]')

ecoli_codon_dict = {
	'A': ['GCG', 'GCC', 'GCA', 'GCT'],
	'C': ['TGC', 'TGT'],
	'D': ['GAT', 'GAC'],
	'E': ['GAA', 'GAG'],
	'F': ['TTT', 'TTC'],
	'G': ['GGC', 'GGT', 'GGG', 'GGA'],
	'H': ['CAT', 'CAC'],
	'I': ['ATT', 'ATC', 'ATA'],
	'K': ['AAA', 'AAG'],
	'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'],
	'M': ['ATG'],
	'N': ['AAC', 'AAT'],
	'P': ['CCG', 'CCA', 'CCT', 'CCC'],
	'Q': ['CAG', 'CAA'],
	'R': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'], 
	'S': ['AGC', 'TCT', 'AGT', 'TCC', 'TCA', 'TCG'],
	'T': ['ACC', 'ACG', 'ACT', 'ACA'],
	'V': ['GTG', 'GTT', 'GTC', 'GTA'],
	'W': ['TGG'],
	'Y': ['TAT', 'TAC'],
	'*': ['TAA', 'TGA', 'TAG']
}

saccharomyces_codon_dict = {
	'A': ['GCT', 'GCC', 'GCA', 'GCG'],
	'C': ['TGT', 'TGC'],
	'D': ['GAT', 'GAC'],
	'E': ['GAA', 'GAG'],
	'F': ['TTT', 'TTC'],
	'G': ['GGT', 'GGA', 'GGC', 'GGG'],
	'H': ['CAT', 'CAC'],
	'I': ['ATT', 'ATC', 'ATA'],
	'K': ['AAA', 'AAG'],
	'L': ['TTG', 'TTA', 'CTA', 'CTT', 'CTG', 'CTC'],
	'M': ['ATG'],
	'N': ['AAT', 'AAC'],
	'P': ['CCT', 'CCC', 'CCA', 'CCG'],
	'Q': ['CAA', 'CAG'],
	'R': ['AGA', 'AGG', 'CGT', 'CGA', 'CGC', 'CGG'], 
	'S': ['TCT', 'AGT', 'TCA', 'TCC', 'AGC', 'TCG'],
	'T': ['ACT', 'ACA', 'ACC', 'ACG'],
	'V': ['GTT', 'GTA', 'GTC', 'GTG'],
	'W': ['TGG'],
	'Y': ['TAT', 'TAC'],
	'*': ['TAA', 'TAG', 'TGA']
}

pichia_codon_dict = {
	'A': ['GCT', 'GCC', 'GCA', 'GCG'],
	'C': ['TGT', 'TGC'],
	'D': ['GAT', 'GAC'],
	'E': ['GAA', 'GAG'],
	'F': ['TTT', 'TTC'],
	'G': ['GGT', 'GGA', 'GGC', 'GGG'],
	'H': ['CAT', 'CAC'],
	'I': ['ATT', 'ATC', 'ATA'],
	'K': ['AAG', 'AAA'],
	'L': ['TTG', 'CTT', 'CTG', 'TTA', 'CTA', 'CTC'],
	'M': ['ATG'],
	'N': ['AAT', 'AAC'],
	'P': ['CCA', 'CCT', 'CCC', 'CCG'],
	'Q': ['CAA', 'CAG'],
	'R': ['AGA', 'AGG', 'CGT', 'CGA', 'CGC', 'CGG'], 
	'S': ['TCT', 'TCC', 'TCA', 'AGT', 'AGC', 'TCG'],
	'T': ['ACT', 'ACA', 'ACC', 'ACG'],
	'V': ['GTT', 'GTC', 'GTG', 'GTA'],
	'W': ['TGG'],
	'Y': ['TAC', 'TAT'],
	'*': ['TAA', 'TAG', 'TGA']
}

ecoli_codon_freq = {'TTT': 28.9, 'TTC': 18.8, 'TTA': 17.5, 'TTG': 18.6, 'CTT': 12.7, 'CTC': 14.1, 'CTA': 3.4, 'CTG': 54.9, 'ATT': 33.9, 'ATC': 31.0, 'ATA': 5.0, 'ATG': 37.4, 'GTT': 19.6, 'GTC': 14.3, 'GTA': 10.6, 'GTG': 33.9, 'TCT': 8.5, 'TCC': 8.0, 'TCA': 6.1, 'TCG': 11.4, 'CCT': 5.8, 'CCC': 2.4, 'CCA': 7.4, 'CCG': 24.9, 'ACT': 7.7, 'ACC': 25.2, 'ACA': 6.1, 'ACG': 14.6, 'GCT': 13.8, 'GCC': 25.5, 'GCA': 19.6, 'GCG': 32.6, 'TAT': 18.6, 'TAC': 8.5, 'TAA': 1.9, 'TAG': 0.3, 'CAT': 9.3, 'CAC': 7.2, 'CAA': 13.5, 'CAG': 24.7, 'AAT': 21.2, 'AAC': 15.9, 'AAA': 29.2, 'AAG': 8.8, 'GAT': 30.0, 'GAC': 15.1, 'GAA': 29.4, 'GAG': 18.0, 'TGT': 4.2, 'TGC': 5.8, 'TGA': 0.8, 'TGG': 12.7, 'CGT': 16.4, 'CGC': 18.8, 'CGA': 2.4, 'CGG': 5.0, 'AGT': 9.0, 'AGC': 14.3, 'AGA': 2.4, 'AGG': 2.1, 'GGT': 24.4, 'GGC': 33.1, 'GGA': 8.2, 'GGG': 14.3} 
#codon_freq amounts from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=413997

genscript_codon_freq = {'TTT': 22.1, 'TTC': 16.0, 'TTA': 14.3, 'TTG': 13.0, 'CTT': 11.9, 'CTC': 10.2, 'CTA': 4.2, 'CTG': 48.4, 'ATT': 29.8, 'ATC': 23.7, 'ATA': 6.8, 'ATG': 26.4, 'GTT': 19.8, 'GTC': 14.3, 'GTA': 11.6, 'GTG': 24.4, 'TCT': 10.4, 'TCC': 9.1, 'TCA': 8.9, 'TCG': 8.5, 'CCT': 7.5, 'CCC': 5.4, 'CCA': 8.6, 'CCG': 20.9, 'ACT': 10.3, 'ACC': 22.0, 'ACA': 9.3, 'ACG': 13.7, 'GCT': 17.1, 'GCC': 24.2, 'GCA': 21.2, 'GCG': 30.1, 'TAT': 17.5, 'TAC': 12.2, 'TAA': 2.0, 'TAG': 0.3, 'CAT': 12.5, 'CAC': 9.3, 'CAA': 14.6, 'CAG': 28.4, 'AAT': 20.6, 'AAC': 21.4, 'AAA': 35.3, 'AAG': 12.4, 'GAT': 32.7, 'GAC': 19.2, 'GAA': 39.1, 'GAG': 18.7, 'TGT': 5.2, 'TGC': 6.1, 'TGA': 1.0, 'TGG': 13.9, 'CGT': 20.0, 'CGC': 19.7, 'CGA': 3.8, 'CGG': 5.9, 'AGT': 9.9, 'AGC': 15.2, 'AGA': 3.6, 'AGG': 2.1, 'GGT': 25.5, 'GGC': 27.1, 'GGA': 9.5, 'GGG': 11.3} 
#genscript_codon_freq amounts from https://www.genscript.com/tools/codon-frequency-table, can use by editing codon_frequency

saccharomyces_codon_freq = {'TTT': 26.1, 'TTC': 18.2, 'TTA': 26.4, 'TTG': 27.1, 'TAT': 18.8, 'TAC': 14.7, 'TAA': 1.0, 'TAG': 0.5, 'CTT': 12.2, 'CTC': 5.4, 'CTA': 13.4, 'CTG': 10.4, 'CAT': 13.7, 'CAC': 7.8, 'CAA': 27.5, 'CAG': 12.2, 'ATT': 30.2, 'ATC': 17.1, 'ATA': 17.8, 'ATG': 20.9, 'AAT': 36.0, 'AAC': 24.9, 'AAA': 42.2, 'AAG': 30.8, 'GTT': 22.0, 'GTC': 11.6, 'GTA': 11.8, 'GTG': 10.6, 'GAT': 37.8, 'GAC': 20.3, 'GAA': 45.9, 'GAG': 19.1, 'TCT': 23.6, 'TCC': 14.2, 'TCA': 18.8, 'TCG': 8.6, 'TGT': 8.0, 'TGC': 4.7, 'TGA': 0.6, 'TGG': 1.00, 'CCT': 13.6, 'CCC': 6.8, 'CCA': 18.2, 'CCG': 5.3, 'CGT': 6.5, 'CGC': 2.6, 'CGA': 3.0, 'CGG': 1.7, 'ACT': 20.3, 'ACC': 12.6, 'ACA': 17.8, 'ACG': 7.9, 'AGT': 14.2, 'AGC': 9.7, 'AGA': 21.3, 'AGG': 9.2, 'GCT': 21.1, 'GCC': 12.5, 'GCA': 16.2, 'GCG': 6.1, 'GGT': 23.9, 'GGC': 9.7, 'GGA': 10.9, 'GGG': 6.0}
#https://www.genscript.com/tools/codon-frequency-table

pichia_codon_freq = {
	'TTT': 23.9, 'TTC': 19.1, 'TTA': 14.9, 'TTG': 31.4, 
	'TAT': 14.7, 'TAC': 18.3, 'TAA': 0.9, 'TAG': 0.5, 
	'CTT': 16.0, 'CTC': 7.6, 'CTA': 11.2, 'CTG': 15.3, 
	'CAT': 10.5, 'CAC': 8.9, 'CAA': 23.9, 'CAG': 14.5, 
	'ATT': 31.7, 'ATC': 19.3, 'ATA': 11.5, 'ATG': 19.2, 
	'AAT': 23.5, 'AAC': 25.7, 'AAA': 30.2, 'AAG': 34.4, 
	'GTT': 26.7, 'GTC': 14.5, 'GTA': 10.1, 'GTG': 12.8, 
	'GAT': 37.2, 'GAC': 26.2, 'GAA': 40.2, 'GAG': 29.6, 
	'TCT': 23.5, 'TCC': 16.3, 'TCA': 15.6, 'TCG': 7.2, 
	'TGT': 8.3, 'TGC': 4.5, 'TGA': 0.3, 'TGG': 9.9, 
	'CCT': 15.3, 'CCC': 6.7, 'CCA': 17.1, 'CCG': 4.1, 
	'CGT': 6.9, 'CGC': 2.3, 'CGA': 4.6, 'CGG': 2.2, 
	'ACT': 23.3, 'ACC': 13.7, 'ACA': 14.3, 'ACG': 6.3, 
	'AGT': 12.1, 'AGC': 7.4, 'AGA': 19.9, 'AGG': 6.6, 
	'GCT': 29.6, 'GCC': 16.7, 'GCA': 15.9, 'GCG': 3.7, 
	'GGT': 26.6, 'GGC': 8.6, 'GGA': 20.0, 'GGG': 6.4}


if args.species == 'pichia':
	codon_dict = pichia_codon_dict
	codon_freq = pichia_codon_freq
elif args.species == 'ecoli':
	codon_dict = pichia_codon_dict
	codon_freq = genscript_codon_freq
else:
	sys.exit('Unknown species!')

def max_frequency(codon):
	for aa in codon_dict:
			if codon in codon_dict[aa]:
				all_codons = codon_dict[aa]
				max_score = 0.0
				for codon in all_codons:
					if codon_freq[codon] >= max_score:
						max_score = codon_freq[codon]
				return max_score


def calculate_cai(codon_sequence):
	codon_list = [codon_sequence[i:i+3] for i in range(0, len(codon_sequence), 3)]
	codon_score = 1.0
	for codon in codon_list:
		max_score = max_frequency(codon)
		codon_score *= (codon_freq[codon] / max_score)
	geom_mean = codon_score ** (1. / float(len(codon_list)))
	return(geom_mean)


for file in fasta_files:
	with open(file, 'r') as f:
		# output_file = open(str(sys.argv[1]).replace('.fasta', '_codons.fasta'), 'w')
		lines = f.readlines()
		for index, line in enumerate(lines):
			if line[0] == '>': 
				new_index = index + 1
				try:
					if lines[new_index][0] == '>':
						sys.exit('Fasta file must be in format:\n >name\n amino_acid seq')
				except:
					sys.exit('Fasta file must be in format:\n >name\n amino_acid seq')
				codon_sequence = ''
				try:
					while lines[new_index][0] != '>':
						while lines[new_index][0] == '#':
							new_index += 1
						codon_sequence += lines[new_index].replace('\n', '')
						new_index += 1
						if lines[new_index][0] == '>':
							pass
				except:
					pass

				print(line.replace('\n', ''))
				cai = calculate_cai(codon_sequence)
				print(round(cai, 2))


