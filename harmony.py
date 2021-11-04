from rt_data import get_codon_dict
import sys
import argparse

parser = argparse.ArgumentParser(description='Manage script input')
parser.add_argument("fasta_files", default=None, nargs="*", help="Path to the input fasta file(s).")
parser.add_argument("-s", "--species", default='ecoli', help="Expression species to optimize for. Currently only accepts ecoli and pichia")
parser.add_argument("-o", "--opt", default='cai', help="Do you want to maximize CAI, or achieve codon harmony? Currently only accepts cai (default) or harmony")
# parser.add_argument("-o", "--outpath", default='.', help="Path to output folder")
args = parser.parse_args()

if len(args.fasta_files) >= 1:
    fasta_files = args.fasta_files
else:
    # files = glob.glob('.fasta')
    sys.exit('You have to provide a fasta file to optimize!\nUsage: python reverse_translate.py FASTA_FILE [optional: --species pichia]')

# get the correct codon usage tables
try:
	codon_dict, codon_freq = get_codon_dict(args.species)
except:
	sys.exit('Unknown species!')

opt = args.opt

def get_aa_from_codon(codon):
	#returns the aa corresponding to the input codon
	for aa in codon_dict:
		if codon in codon_dict[aa]:
			return aa

def translate(dna_seq):

	translation = []
	lst = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
	for codon in lst:
		translation.append(get_aa_from_codon(codon))
	return ''.join(translation)

def calculate_harmony(codon_sequence):

	translated_seq = translate(codon_sequence)
	
	# calculate "codon harmony", i.e. how closely the current sequence matches the distribution of host genome codons
	# penalty is calculated via summing the deviations of the usage in the optimized sequence vs. the host genome
	# for 20 amino acids, the sum of codon usage corresponding to each amino acid is 1, so the total penalty is max. 20
	# to make this more readible, let's multiply the penalty by 5 (base 100), delete from 100, and express as a %

	harmony_score = 0
	codon_list = [codon_sequence[i:i+3] for i in range(0, len(codon_sequence), 3)]

	all_codons = codon_freq.keys()


	for c in all_codons:
		# get the sum of all codon frequency scores for this amino acid
		current_aa = get_aa_from_codon(c)
		
		if translated_seq.count(current_aa) != 0:
			all_codons_from_current_aa = codon_dict[current_aa]
			total_codon_freq = float(sum(codon_freq[i] for i in all_codons_from_current_aa))
			
			# get the target frequency
			target_freq = codon_freq[c] / total_codon_freq

			# get the actual frequency
			actual_num_codons = codon_list.count(c)
			actual_total_codons_from_aa = translated_seq.count(current_aa)
			actual_freq = float(actual_num_codons) / float(actual_total_codons_from_aa)

			difference = abs(actual_freq - target_freq)
			harmony_score += difference

			print(c, actual_num_codons, actual_freq, target_freq)

	return(100 - (harmony_score * 5.))

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
				codon_seq = ''
				try:
					while lines[new_index][0] != '>':
						while lines[new_index][0] == '#':
							new_index += 1
						codon_seq += lines[new_index].replace('\n', '')
						new_index += 1
						if lines[new_index][0] == '>':
							pass
				except:
					pass

				print(line.replace('\n', ''))
				print(codon_seq)
				harmony = calculate_harmony(codon_seq)
				print(round(harmony, 2))