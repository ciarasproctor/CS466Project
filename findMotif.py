#!/usr/bin/python

import argparse, random, numpy, itertools, os
from pyfasta import Fasta

_debugging_endabled_ = False

def debug(string):
	if _debugging_endabled_:
		print string

def weighted_choice(weights):
	totals = []
	running_total = 0

	for w in weights:
		running_total += w
		totals.append(running_total)

		rnd = random.random() * running_total
		for i, total in enumerate(totals):
			if rnd < total:
				return i

class NGram:
	def __init__(self, sequence, n):
		self.position = iter(sequence)
		self.size = n

	def __iter__(self):
		return self

	def next(self):
		ngram = []
		for c in itertools.islice(self.position, 0, self.size):
			ngram.append(c)

		if len(acid) != self.size:
			raise StopIteration

		return ngram

class MotifFinder:
	def __init__(self, sequences, motif_length):
		self.bases = 4
		self.base_map = { 'g': 0, 'c': 1, 'a': 2, 't': 3 }
		self.pseudocounts = [ .5, .5, .5, .5 ]
		self.pseudocounts_total = 0
		self.motif_length = motif_length
		self.sequences = sequences
		self.sequence_count = len(sequences.keys())
		self.residue_frequencies = numpy.zeros((self.motif_length, self.bases))
		self.background_frequencies = numpy.zeros((self.bases))
		self.alignments = numpy.zeros(self.sequence_count)

		for c in self.pseudocounts:
			self.pseudocounts_total += c

		# initialize the alignments randomly
		for (i, a) in enumerate(self.alignments):
			self.alignments[i] = random.randint(0, len(self.sequence(i)) - 1)

	def sequence(self, i):
		key = self.sequences.keys()[i]
		return self.sequences[key]

	def base_index(self, base):
		return self.base_map[base.lower()]

	def gibbs_sampler(self):
		z = random.randint(0, self.sequence_count - 1)
		sequence = self.sequence(z)
		self.predictive_update(z)

		self.alignments[z] = self.sample(self.sequence(z))

		debug(self)

	def predictive_update(self, z):
		counts = numpy.zeros((self.motif_length, self.bases))
		background_counts = numpy.zeros((self.bases))
		background_total = 0

		for s in range(0, self.sequence_count):
			if s != z:
				start = self.alignments[s]
				sequence = self.sequence(s)
				end = self.motif_length + start - 1
				for (i, base) in enumerate(sequence):
					j = self.base_index(base)
					if i < start or i > end:
						background_counts[j] += 1
						background_total += 1
					else:
						counts[i - start, j] += 1

		# update the residue frequencies
		for i in range(0, self.motif_length):
			for j in range(0, self.bases):
				self.residue_frequencies[i,j] = (counts[i,j] + self.pseudocounts[j]) / (self.sequence_count - 1 + self.pseudocounts_total)

		background_length = 0
		for (i, (k, seq)) in enumerate(self.sequences.iteritems()):
			if i != z:
				background_length += len(seq) - self.motif_length

		# update the background frequencies
		for i in range(0, self.bases):
			self.background_frequencies[i] = (background_counts[i] + self.pseudocounts[i]) / (background_length + self.pseudocounts_total)


	def sample(self, sequence):
		max_start = len(sequence) - self.motif_length
		likelihoods = numpy.zeros((max_start))
		likelihoods_total = 0
		for start in range(0, max_start):
			likelihoods[start] = self.likelihood(sequence, start)
			likelihoods_total += likelihoods[start]

		# normalize the likelihoods
		multiplier = 1 / likelihoods_total
		for likelihood in likelihoods:
			likelihood *= multiplier

		return weighted_choice(likelihoods)

	def likelihood(self, sequence, start):
		motif_probability = 1
		background_probability = 1
		for (i, base) in enumerate(itertools.islice(sequence, start, self.motif_length)):
			base_index = self.base_index(base)
			motif_probability *= self.residue_frequencies[i, base_index]
			background_probability *= self.background_frequencies[base_index]

		return motif_probability / background_probability

	def __str__(self):
		sequence_strs = []
		for (k, seq) in self.sequences.iteritems():
			sequence_strs.append(str(seq))
		str_sequences = '\"'
		str_sequences += ('\"\n\t%s\"' % (' ' * len('sequences[%d] = ' % self.sequence_count))).join(sequence_strs)
		str_sequences += '\"'
		residue_strs = []
		for line in self.residue_frequencies:
			residue_strs.append(str(line))
		str_residue = '['
		str_residue += ('\n\t%s' % (' ' * len('residue_frequencies = ['))).join(residue_strs)
		str_residue += ']'
		return "MotifFinder {\n\tmotif_length = %s\n\tsequences[%d] = %s\n\tresidue_frequencies = %s\n\tbackground_frequences = %s\n\talignments = %s\n}" % (self.motif_length, self.sequence_count, str_sequences, str_residue, self.background_frequencies, self.alignments)

	def __repr__(self):
		return self.__str__()

def main():
	global _debugging_endabled_
	parser = argparse.ArgumentParser(description='Motif finding using the Gibbs Motif Sampling Algorithm', formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='''Generates the following files (in directory):
  * predictedmotif.txt
  * predictedsites.txt''')
	parser.add_argument('directory', type=str, help='Read and write files to this directory', default='.')
	parser.add_argument('iterations', type=int, help='Number of Gibbs Sampler iterations to perform')
	parser.add_argument('--debug', action='store_true', help='Turn on debugging information', default=False)
	args = parser.parse_args()

	_debugging_endabled_ = args.debug
	length_path = os.path.join(args.directory, 'motiflength.txt')
	sequences_path = os.path.join(args.directory, 'sequences.fa')

	debug("using length from: %s" % length_path)
	debug("using sequences from: %s" % sequences_path)

	length_file = open(length_path, 'r')
	motif_length = int(length_file.readline())

	sequences = Fasta(sequences_path)

	finder = MotifFinder(sequences, motif_length)
	debug(finder)

	for i in range(0, args.iterations):
		finder.gibbs_sampler()

	print "Motif: %s" % finder.alignments

if __name__ == '__main__':
	main()
