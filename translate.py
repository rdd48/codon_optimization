import sys

file = sys.argv[1]

codon_dict = {
#AGG, AGA, CGG, CGA (R), GGA (G), ATA (I), CTA (L), CCC (P), ACA (T) associated with poor expression in E. coli. Indexed by frequency found in E. coli B
	'A': ['GCG', 'GCC', 'GCA', 'GCT'],
	'C': ['TGC', 'TGT'],
	'D': ['GAT', 'GAC'],
	'E': ['GAA', 'GAG'],
	'F': ['TTT', 'TTC'],
	'G': ['GGC', 'GGT', 'GGG', 'GGA'],
	'H': ['CAT', 'CAC'],
	'I': ['ATT', 'ATC', 'ATA'],
	'K': ['AAA', 'AAG'],
	'L': ['CTG', 'TTG', 'TTA', 'CTC', 'CTT', 'CTA'],
	'M': ['ATG'],
	'N': ['AAT', 'AAC'],
	'P': ['CCG', 'CCA', 'CCT', 'CCC'],
	'Q': ['CAG', 'CAA'],
	'R': ['CGC', 'CGT', 'AGG', 'AGA', 'CGG', 'CGA'], 
	'S': ['AGC', 'TCG', 'AGT', 'TCT', 'TCC', 'TCA'],
	'T': ['ACC', 'ACG', 'ACT', 'ACA'],
	'V': ['GTG', 'GTT', 'GTC', 'GTA'],
	'W': ['TGG'],
	'Y': ['TAT', 'TAC'],
	'*': ['TAA', 'TAG', 'TGA']
}

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

with open(file, 'r') as f:
	output_file = open(str(file).replace('.fasta', '_translation.fasta'), 'w')
	lines = f.readlines()
	for index, line in enumerate(lines):
		if line[0] == '>': 
			new_index = index + 1
			try:
				if lines[new_index][0] == '>':
					sys.exit('Fasta file must be in format:\n >name\n dna seq')
			except:
				sys.exit('Fasta file must be in format:\n >name\n dna seq')
			seq_to_translate = ''
			try:
				while lines[new_index][0] != '>':
					while lines[new_index][0] in ['#', '/']:
						new_index += 1
					seq_to_translate += lines[new_index].replace('\n', '')
					new_index += 1
					if lines[new_index][0] == '>':
						pass
			except:
				pass
			seq_to_translate = seq_to_translate.upper()
			output_file.write(str(line))

			final_seq = translate(seq_to_translate)
			output_file.write(final_seq)
			output_file.write('\n')

	output_file.close()
