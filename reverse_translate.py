from rt_data import get_codon_dict
import random
import sys
import argparse
import glob

parser = argparse.ArgumentParser(description='Manage script input')
parser.add_argument("fasta_files", default=None, nargs="*", help="Path to the input fasta file(s).")

parser.add_argument("-s", "--species", default='ecoli', help="Expression species to optimize for. Currently only accepts ecoli and pichia")
parser.add_argument("-o", "--opt", default='cai', help="Do you want to maximize CAI, or achieve codon harmony? Currently only accepts cai (default) or harmony")

# num sequences to print
parser.add_argument("-n", "--num_outputs", default=1, help="How many unique sequences do you want to try for?")

args = parser.parse_args()


if len(args.fasta_files) >= 1:
    fasta_files = args.fasta_files
else:
    all_files = glob.glob('*.fasta')
    fasta_files = [f if '_codons.fasta' not in f else f for f in all_files]
    sys.exit('You have to provide a fasta file to optimize!\nUsage: python reverse_translate.py FASTA_FILE [optional: --species pichia]')

opt = args.opt

# get the correct codon usage tables
try:
	codon_dict, codon_freq = get_codon_dict(args.species)
except:
	sys.exit('Unknown species!')


restriction_enzymes = {
	'bamhi':  'GGATCC',
	'ndei': 'CATATG', 
	'pvuii': 'CAGCTG',
	'saci': 'GAGCTC',
	'xhoi': 'CTCGAG'
}

def get_indices_of_aa_from_seq(aa, sequence):
	indices = []
	for i in range(len(sequence)):
		if sequence[i] == aa:
			indices.append(i)
	return(indices)

def get_codon_and_update_usage_list(codon_usage_list, all_codons):
	for i in range(len(codon_usage_list)):
		if codon_usage_list[i] > 0:
			codon_usage_list[i] -= 1
			return(all_codons[i], codon_usage_list)

def get_harmonic_sequence(sequence):
	'''sequence is a string of amino acids; returns a list of codons.
	Goal: return a codon sequence that "harmonizes" to codon usage in the host genome.
	I can think of only the following, potentially inefficient way
	1. loop through all amino acids (aa)
	2. get all indices in our input sequence for the current aa
	3. get the target number of each codon to use
	   (i.e., if two codons are used 70% and 30% of the time, and there are 10 of those aas,
		  use 7 of the first codon and 3 of the second)
	4. Randomly shuffle these codons throughout the sequence
	5. TBD: remove rare codons from being considered?'''

	aas = "ACDEFGHIKLMNPQRSTVWY*"
	
	for i in sequence:
		if i not in aas:
			exit('Error: sequence contains non-canonical amino acids!')
	# initialize codon sequence with list of 0s of same length:
	codon_seq = [0 for i in range(len(sequence))]

	for aa in aas:
		indices = get_indices_of_aa_from_seq(aa, sequence)
		num_aas = len(indices)

		all_codons_from_current_aa = codon_dict[aa]
		aa_codon_freq = [codon_freq[i] for i in all_codons_from_current_aa]
		
		# get the target frequency
		normalized_codon_freq = [(codon_freq[c] / sum(aa_codon_freq)) for c in all_codons_from_current_aa]

		# get the number of each possible codon to use, e.g.: use x number of GCAs
		# indices of this list correspond to indices of normalized_codon_freq
		codon_nums_to_use = [round(freq * float(num_aas)) for freq in normalized_codon_freq]

		random.shuffle(indices)

		for i in indices:
			# if the round function creates a list that doesn't sum to the total codons,
			# start putting in the most used codon (at the 0 index of the all_codons list)
			if sum(codon_nums_to_use) == 0:
				codon_seq[i] = all_codons_from_current_aa[0]
			else:
				codon, codon_nums_to_use  = get_codon_and_update_usage_list(codon_nums_to_use, all_codons_from_current_aa)
				codon_seq[i] = codon


	return codon_seq 


def first_pass(sequence):
	'''creates the first reverse translation by trying to get perfect codon "harmony", i.e. synchronicity with the host genome.'''
	rev_trans = []
	for aa in sequence:
		if aa not in 'ACDEFGHIKLMNPQRSTVWY*':
			sys.exit('Error: sequence contains non-canonical amino acids!')
		if args.opt == 'cai':
			rev_trans.append(codon_dict[aa][0])
		elif args.opt == 'harmony':
			return get_harmonic_sequence(sequence)
		else:
			sys.exit('Codon optimization method not found! Try: --opt cai or --opt harmony')
	return rev_trans

def create_sliding_windows(dna_seq, window_len):
	'''takes codon list and returns list of n-numbered nucleotide sliding windows'''
	joined_list = ''.join(dna_seq)
	sliding_windows = []
	count = 0
	for nt in joined_list:
		if len(joined_list[count : count + window_len]) == window_len:
			sliding_windows.append(joined_list[count : count + window_len])
			count += 1
	return sliding_windows

def get_codon_pos(nt_pos):
	'''takes in a nucleotide position, outputs the codon position'''
	return int((nt_pos - 1) / 3) + 1

def get_aa_from_codon(codon):
	'''returns the aa corresponding to the input codon'''
	for aa in codon_dict:
		if codon in codon_dict[aa]:
			return aa

def return_new_codon(starting_codon):
	'''generates the next most frequent codon some percent of the time (currently set at 50%),
	otherwise generates the previous codon'''
	if starting_codon == 'M':
		return starting_codon, 'M'
	elif starting_codon == 'W':
		return starting_codon, 'W'
	aa = get_aa_from_codon(starting_codon)
	codon_options = codon_dict[aa]
	codon_index = codon_options.index(starting_codon)
	random_int = random.randint(1,2)
	if random_int == 2:
		try:
			new_codon = codon_options[codon_index + 1]
		except:
			new_codon = codon_options[0]
	else:
		if codon_index != 0:
			new_codon = codon_options[codon_index - 1]
		else:
			new_codon = codon_options[codon_index + 1]
	return new_codon, aa

def replace_codon(starting_codon, low_gc=False, high_gc=False):
	'''replaces codon with the next-most codon chosen by return_new_codon'''
	if starting_codon == 'ATG' or starting_codon == 'TGG':
		return starting_codon
	new_codon, aa = return_new_codon(starting_codon)
	count = 0
	if low_gc == True:
		starting_codon_gc = (starting_codon.count('G') + starting_codon.count('C'))
		new_codon_gc = (new_codon.count('G') + new_codon.count('C'))
		if new_codon_gc < starting_codon_gc:
			return new_codon
		else:
			while new_codon_gc >= starting_codon_gc and count < 5:
				new_codon, aa = return_new_codon(starting_codon)
				new_codon_gc = (new_codon.count('G') + new_codon.count('C'))
				count += 1
			return new_codon
	if high_gc == True:
		starting_codon_gc = (starting_codon.count('G') + starting_codon.count('C'))
		new_codon_gc = (new_codon.count('G') + new_codon.count('C'))
		if new_codon_gc > starting_codon_gc:
			return new_codon
		else:
			while new_codon_gc <= starting_codon_gc and count < 5:
				new_codon, aa = return_new_codon(starting_codon)
				new_codon_gc = (new_codon.count('G') + new_codon.count('C'))
				count += 1
			return new_codon
	return new_codon

def find_duplicates(lst):
	'''returns the duplicate sequence with the window position'''
	first_time = []
	duplicates = []
	count = 1
	for item in lst:
		if item not in first_time:
			first_time.append(item)
			count += 1
		else:
			duplicates.append([item, count])
			count +=1
	return duplicates

def replace_duplicates(codon_sequence, duplicate_len = 8):
	'''takes in a list of codons and outputs a new sequence (with any remaining duplicate positions) 
	where the duplicate positions have been randomly changed'''
	sliding_windows = create_sliding_windows(codon_sequence, duplicate_len)
	duplicates = find_duplicates(sliding_windows)
	for i in duplicates:
		duplicate_codon_pos = get_codon_pos(i[1])
		old_codon = codon_sequence[duplicate_codon_pos]
		new_codon = replace_codon(old_codon)
		codon_sequence[duplicate_codon_pos] = new_codon
	new_duplicates = find_duplicates(create_sliding_windows(codon_sequence, duplicate_len))
	return codon_sequence, new_duplicates

def count_gc(codon_sequence):
	g_count = 0
	c_count = 0
	for codon in codon_sequence:
		for nt in codon:
			if nt == 'G':
				g_count += 1
			elif nt == 'C':
				c_count += 1
	seq_len = float(len(codon_sequence) * 3)
	return(float(g_count + c_count) / seq_len)

def find_high_gc(sliding_windows, gc_cutoff = 0.65):
	'''returns the windows with %gc above a specified cutoff'''
	count = 0
	high_gc_windows = []
	best_gc_window = []
	for window in sliding_windows:
		gc_in_window = window.count('G') + window.count('C')
		gc_ratio = float(gc_in_window) / float(len(window))
		count += 1
		if gc_ratio >= gc_cutoff:
			high_gc_windows.append([window, count])
	return high_gc_windows

def find_low_gc(sliding_windows, gc_cutoff = 0.35):
	'''returns the windows with %gc below a specified cutoff'''
	count = 0
	low_gc_windows = []
	for window in sliding_windows:
		gc_in_window = window.count('G') + window.count('C')
		gc_ratio = float(gc_in_window) / float(len(window))
		count += 1
		if gc_ratio <= gc_cutoff:
			low_gc_windows.append([window, count])
	return low_gc_windows

def reduce_gc(codon_sequence, window_length, gc_cutoff = 0.65):
	'''attempts to reduce gc content by replacing higher gc codons with lower gc codons'''
	sliding_windows = create_sliding_windows(codon_sequence, window_length)
	high_gc_windows = find_high_gc(sliding_windows, gc_cutoff)
	random_int = random.randint(1,10)
	for i in high_gc_windows[random_int::5]: #reduce gc for every 10th window, starting at a random position between 1-10
		first_codon_pos = get_codon_pos(i[1])
		old = codon_sequence[first_codon_pos - 1]
		new = replace_codon(old, low_gc = True)
		if (old.count('G') + old.count('C')) > (new.count('G') + new.count('C')):
			codon_sequence[first_codon_pos - 1] = new
	return codon_sequence, high_gc_windows

def increase_gc(codon_sequence, window_length, gc_cutoff = 0.30):
	'''attempts to increase gc content by replacing lower gc codons with higher gc codons'''
	sliding_windows = create_sliding_windows(codon_sequence, window_length)
	low_gc_windows = find_low_gc(sliding_windows, gc_cutoff)
	random_int = random.randint(1,10)
	for i in low_gc_windows[random_int::5]: #reduce gc for every 10th window, starting at a random position between 1-10
		first_codon_pos = get_codon_pos(i[1])
		old = codon_sequence[first_codon_pos - 1]
		new = replace_codon(old, high_gc = True)
		if (old.count('G') + old.count('C')) < (new.count('G') + new.count('C')):
			codon_sequence[first_codon_pos - 1] = new
	return codon_sequence, low_gc_windows

def reduce_gc_last_pos(codon_sequence, window_length, gc_cutoff = 0.65):
	'''attempts to reduce gc content by replacing higher gc codons with lower gc codons'''
	sliding_windows = create_sliding_windows(codon_sequence, window_length)
	high_gc_windows = find_high_gc(sliding_windows, gc_cutoff)
	already_mutated = []
	for i in high_gc_windows:
		#first_codon_pos = get_codon_pos(i[1]) - 1
		last_codon_pos = get_codon_pos(i[1] + window_length - 1) - 1
		old = codon_sequence[last_codon_pos]
		if last_codon_pos not in already_mutated:
			already_mutated.append(last_codon_pos)
			new = replace_codon(old, low_gc = True)
			if (old.count('G') + old.count('C')) > (new.count('G') + new.count('C')):
				codon_sequence[last_codon_pos] = new
	return codon_sequence, high_gc_windows

def codon_frequency(codon_sequence):
	'''returns a 'score' based off how often that codon is found in e. coli'''
	total_score = 0
	for codon in codon_sequence:
		total_score += codon_freq[codon]
	return total_score

def find_best_expressing_seq(sequences):
	'''determines best codon_frequency-scoring sequence. currently only really works if all sequences are the same length, which is fine for now.'''
	best_seq = sequences[0]
	best_score = 0

	count = 0

	for seq in sequences:
		current_score = codon_frequency(sequences[count])
		if current_score > best_score:
			best_score = current_score
			best_seq = sequences[count]
		count += 1
	return best_seq

def change_restriction_sites(codon_sequence, enzymes):
	windows = create_sliding_windows(codon_sequence, 6)
	codon_pos = []
	for enzyme in enzymes:
		enzyme = enzyme.lower()
		if enzyme not in restriction_enzymes:
			continue
		else:
			restriction_seq = restriction_enzymes[enzyme]
			count = 0
			for window in windows:
				count += 1
				if window == restriction_seq:
					codon_pos.append(get_codon_pos(count))
			if codon_pos != []:
				for pos in codon_pos:
					new_codon = replace_codon(codon_sequence[pos - 1])
					codon_sequence[pos - 1] = new_codon
	return codon_sequence

def replace_negative_cis_element(seq_list):
	seq_str = ''.join(seq_list)
	splice_sites = ['GGTAAG', 'GGTGAT', 'ATTTTTTA', 'AAAAAAA']
	
	for ss in splice_sites:
		for i in range(len(seq_str) - len(ss) + 1):
			if seq_str[i:i+len(ss)] == ss:
				codon_idx = get_codon_pos(i + 1) - 1 # workaround, but this should've been indexed by zero when i first wrote it...
				seq_list[codon_idx] = return_new_codon(seq_list[codon_idx])[0]
				
		seq_str = ''.join(seq_list)
	
	return [seq_str[i:i+3] for i in range(0, len(seq_str), 3)]

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

def calculate_harmony(codon_sequence):
	'''calculate "codon harmony", i.e. how closely the current sequence matches the distribution of host genome codons
	penalty is calculated via summing the deviations of the usage in the optimized sequence vs. the host genome
	for 20 amino acids, the sum of codon usage corresponding to each amino acid is 1, so the total penalty is max. 20
	to make this more readible, let's multiply the penalty by 5 (base 100), delete from 100, and express as a %'''

	harmony_score = 0

	codon_list = [codon_sequence[i:i+3] for i in range(0, len(codon_sequence), 3)]

	all_codons = codon_freq.keys()


	for c in all_codons:
		# get the sum of all codon frequency scores for this amino acid
		current_aa = get_aa_from_codon(c)
		
		if seq_to_translate.count(current_aa) != 0:
			all_codons_from_current_aa = codon_dict[current_aa]
			total_codon_freq = float(sum([codon_freq[i] for i in all_codons_from_current_aa]))
			
			# get the target frequency
			target_freq = codon_freq[c] / total_codon_freq

			# get the actual frequency
			actual_num_codons = codon_list.count(c)
			actual_total_codons_from_aa = seq_to_translate.count(current_aa)
			actual_freq = float(actual_num_codons) / float(actual_total_codons_from_aa)

			difference = abs(actual_freq - target_freq)
			harmony_score += difference
		
	return(100 - (harmony_score * 5.))

def find_most_harmonious_seq(sequences, aa_seq):
    """determines best harmony-scoring sequence.
    currently only really works if all sequences are the same length, which is fine for now."""
    best_seq = sequences[0]
    best_harmony = 0

    for seq in sequences:
        current_score = calculate_harmony("".join(seq), aa_seq)
        if current_score > best_harmony:
            best_harmony = current_score
            best_seq = seq
    return best_seq

def optimize_sequence(aa_sequence):
	codon_sequence = first_pass(aa_sequence)
	lower_duplicate_seq, duplicates = replace_duplicates(codon_sequence, 8)
	
	count = 0
	while duplicates != [] and count < 3:
		lower_duplicate_seq, duplicates = replace_duplicates(lower_duplicate_seq, 8)
		#lower_duplicate_seq, duplicates = replace_duplicates(lower_duplicate_seq, 12)
		count += 1
	return lower_duplicate_seq

def main(aa_sequence, num_times_to_loop = 50):
	'''makes multiple (# depends on num_times variable) reverse translations and returns the most codon-optimized. Breaks duplicates one last time after that.'''
	rev_trans_sequences = []
	for i in range(num_times_to_loop):
		optimized_seq = optimize_sequence(aa_sequence)
		
		gc_ratio = count_gc(optimized_seq)
		count = 0

		gc_seq_len = int(len(aa_sequence)/ 10)
		
		if gc_seq_len < 8:
			gc_seq_len = 8

		while count <= 2:
			if gc_ratio <= 0.4:
				optimized_seq, low_gc_windows = increase_gc(optimized_seq, 100, 0.4)
			elif gc_ratio >= 0.6:
				optimized_seq, high_gc_windows = reduce_gc(optimized_seq, 100, 0.6)
			count += 1
		optimized_seq, duplicates = replace_duplicates(optimized_seq, 16)

		optimized_seq = replace_negative_cis_element(optimized_seq)

		rev_trans_sequences.append(optimized_seq)

	optimized_seq = find_best_expressing_seq(rev_trans_sequences)

	#return(''.join(first_pass(aa_sequence)))
	return(''.join(optimized_seq))
	

if __name__ == '__main__':
	for file in fasta_files:
		with open(file, 'r') as f:
			output_file = open(str(file).replace('.fasta', '_codons.fasta'), 'w')
			lines = f.readlines()
			for index, line in enumerate(lines):
				if line[0] == '>': 
					new_index = index + 1
					try:
						if lines[new_index][0] == '>':
							sys.exit('Fasta file must be in format:\n >name\n amino_acid seq')
					except:
						sys.exit('Fasta file must be in format:\n >name\n amino_acid seq')
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

					for i in range(int(args.num_outputs)):
						final_seq = main(seq_to_translate, 20)
						output_file.write(final_seq)
						output_file.write('\n')
						print(line.replace('\n', ''))
						print('CAI: ' + str(calculate_cai(final_seq)))
						print('% GC: ' + str(count_gc([final_seq[i:i+3] for i in range(0, len(final_seq), 3)])))
						print('% Harmony: ' + str(calculate_harmony(final_seq)))

			output_file.close()

