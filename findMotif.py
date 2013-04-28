#!/usr/bin/python

import argparse, random, numpy, itertools
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
		self.pseudocounts = { .5, .5, .5, .5 }
		self.pseudocounts_total = 0
		self.motif_length = motif_length
		self.sequences = sequences
		self.sequence_count = len(sequences.keys())
		self.residue_frequencies = numpy.zeros((self.motif_length, self.bases))
		self.background_frequencies = numpy.zeros(self.bases)
		self.alignments = numpy.zeros(self.sequence_count)

		for c in self.pseudocounts:
			self.pseudocounts_total += c

		# initialize the alignments randomly
		for (i, a) in enumerate(self.alignments):
			a = random.randint(0, len(self.sequences[i]) - 1)

	def gibbs_sampler(self):
		z = random(0, self.sequence_count)
		sequence = self.sequences[z]
		self.predictive_update(z)

		self.alignments[z] = self.sample(z)

		debug(self)

	def predictive_update(self, z):
		counts = numpy.zeros((self.motif_length, self.bases))
		background_counts = numpy.zeros(self.bases)
		background_total = 0
		for s in range(0, self.sequence_count):
			if s != z:
				start = self.alignments[s]
				sequence = self.sequences[s]
				end = self.motif_length + start - 1
				for (i, base) in enumerate(sequence):
					j = base_index(base)
					if i < start or i > end:
						background_counts[j] += 1
						background_total += 1
					else:
						counts[i, j] += 1

		# update the residue frequencies
		for i in range(0, self.motif_length):
			for j in range(0, self.bases):
				self.residue_frequencies[i,j] = (counts[i,j] + self.pseudocounts[j]) / (self.sequence_count - 1 + self.pseudocounts_total)

		background_length = 0
		for (i, seq) in enumerate(self.sequences):
			if i != z:
				background_length += len(seq) - self.motif_length

		# update the background frequencies
		for i in range(0, self.bases):
			self.background_frequencies[i] = (background_counts[i] + self.pseudocounts[i]) / (background_length + self.pseudocounts_total)


	def sample(self, sequence):
		max_start = len(sequence) - self.motif_length
		likelihoods = numpy.zeros(max_start)
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
		motif_probability = 0
		background_probability = 0
		for (i, base) in enumerate(islice(sequence, start, self.motif_length)):
			base_index = base_index(base)
			motif_probability *= self.residue_frequencies[i, base_index]
			background_probability *= self.background_frequencies[base_index]

		return motif_probability / background_probability

	def __str__(self):
		return "MotifFinder {\n\tmotif_length = %s\n\tsequences[%d] = %s\n\tresidue_frequencies = %s\n\tbackground_frequences = %s\n\talignments = %s\n}" % (self.motif_length, self.sequence_count, self.sequences, self.residue_frequencies, self.background_frequencies, self.alignments)

	def __repr__(self):
		return self.__str__()

def main():
	parser = argparse.ArgumentParser(description='Motif finding using the Gibbs Motif Sampling Algorithm', formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='''Generates the following files (optionally prefixed by prefix):
  * predictedmotif.txt
  * predictedsites.txt''')
	parser.add_argument('sequences', type=file, help='The file containing the sequences (FASTA format)')
	parser.add_argument('length', type=file, help='Path to the file containing the motif length')
	parser.add_argument('--prefix', type=str, help='Prefix output files with this string', default='')
	parser.add_argument('--debug', action='store_true', help='Turn on debugging information', default=False)
	args = parser.parse_args()

	_debugging_endabled_ = args.debug
	motif_length = int(args.length.readline())
	sequences = Fasta(args.sequences)

	finder = MotifFinder(sequences, motif_length)
	debug(finder)

	for i in range(0, 10):
		finder.gibbs_sampler()

	print "Motif: %s" % finder

if __name__ == '__main__':
	main()
