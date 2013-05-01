#!/usr/bin/python

import argparse,numpy,os,time,random
from pyfasta import Fasta

def main():
	start_time = time.time()

	parser = argparse.ArgumentParser(description="Motif finding using a conscensus algorithm", formatter_class=argparse.RawDescriptionHelpFormatter, epilog= '''Generates the following files (in directory):
	* predictedmotif.txt
	* predictedsites.txt
	* runningtime.txt''')
	parser.add_argument('directory', type=str, help='Read and write files to this directory', default = '.')
	args = parser.parse_args()

	(sequences,ml) = load_data(args.directory)
	(motif, sites) = find_motif(sequences,ml)
	write_motif(motif,ml,args.directory)
	write_sites(sites,args.directory)
	duration = time.time() - start_time
	write_running_time(duration, args.directory)

def find_motif(sequences,ml):
	seq = ""
	index = 0
	motif = numpy.zeros(shape=(ml,4))
	sites = []
	(motif,loc) = generate_initial_motif(sequences[0],sequences[1],ml,motif)
	sites.append(loc[0])
	sites.append(loc[1])
	for i in range(2,len(sequences)):
		(motif,loc) = find_site(motif,sequences[i],ml)
		sites.append(loc)
	
	return(motif,sites)

def convert_sequences(sequences_d):
	sequences = {}
	for key in sequences_d.keys():
		sequences[int(key)] = sequences_d[key][:]
	return sequences

def generate_initial_motif(seq1,seq2,ml,motif):
	best_motif = []
	best_locations = []
	best_score = 0
	possible_starts = len(seq1) - 1 - ml
	for i in range(possible_starts):
		motif = numpy.zeros(shape=(ml,4))
		l1 = seq1[i:i+ml]
		update_motif(l1,motif)
		for j in range(possible_starts):
			l2 = seq2[j:j+ml]
			cs = consensus_score(l2,motif)
			if cs > best_score:
				best_score = cs
				best_motif = []
				best_locations = []
			if cs == best_score:
				good_motif = numpy.copy(motif)
				update_motif(l2,good_motif)
				best_motif.append(good_motif)
				best_locations.append((i,j))
	random_index = random.randint(0,len(best_motif)-1)
	return (best_motif[random_index], best_locations[random_index])

def update_motif(seq,motif):
	bases = "ACGT"
	for i in range(len(seq)):
		motif[i,bases.index(seq[i])] += 1

def find_site(motif,seq,ml):
	best_score = 0
	best_motifs = []
	best_locations = []
	for i in range(len(seq)-1-ml):
		l = seq[i:i+ml]
		cs = consensus_score(l,motif)
		if cs > best_score:
			best_score = cs
			best_motifs = []
			best_locations = []
		if cs == best_score:
			good_motif = numpy.copy(motif)
			update_motif(l,good_motif)
			best_motifs.append(good_motif)
			best_locations.append(i)
	random_index = random.randint(0,len(best_motifs)-1)
	return (best_motifs[random_index], best_locations[random_index])

def consensus_score(l,motif):
	cs = 0
	bases = "ACGT"
	for i in range(len(l)):
		cs += motif[i,bases.index(l[i])]
	return cs
	
def load_data(directory):
	length_path = os.path.join(directory,"motiflength.txt")
	sequences_path = os.path.join(directory,"sequences.fa")
	length_file = open(length_path,"r")
	ml = int(length_file.readline())
	sequences = []
	sequence_file = open(sequences_path,'r')
	sequence_lines = sequence_file.readlines()
	for line in sequence_lines:
		if len(line) > 0 and line[0] != '>':
			sequences.append(line.strip())
	return(sequences,ml)

def write_motif(pwm,ml,directory):
	f = open(os.path.join(directory,'predictedmotif_b.txt'),'w')
	f.write('>motif %d \n' %ml)
	for row in pwm:
		rowText = ""
		for cell in row:
			rowText += str(cell) + " "
		rowText=rowText.strip()
		f.write(rowText+'\n')
	f.write('<\n')
	f.close()

def write_sites(sites,prefix):
	f = open(os.path.join(prefix,'predictedsites_b.txt'),'w')
	for site in sites:
		f.write(str(site)+'\n')
	f.close()

def write_running_time(duration,directory):
	f = open(os.path.join(directory,"running_time_b.txt"), "w")
	f.write(str(duration) + "\n")
	f.close()

if __name__ == "__main__":
	main()
