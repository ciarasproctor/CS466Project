#!/usr/bin/python

import argparse, random, numpy
from pyfasta import Fasta

class MotifFinder:
	def __init__(self, sequences, motif_length):
		self.motif_length = motif_length
		self.sequences = sequences
		self.sequence_count = len(sequences.keys())
		self.residue_frequencies = numpy.zeros((motif_length, 20))
		self.background_frequenceis = numpy.zeros(20)
		self.alignments = numpy.zeros(self.sequence_count)


def main():
	parser = argparse.ArgumentParser(description='Motif finding using the Gibbs Motif Sampling Algorithm', formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='''Generates the following files (optionally prefixed by prefix):
  * predictedmotif.txt
  * predictedsites.txt''')
	parser.add_argument('sequences', type=file, help='The file containing the sequences (FASTA format)')
	parser.add_argument('length', type=file, help='Path to the file containing the motif length')
	parser.add_argument('--prefix', type=str, help='Prefix output files with this string', default='')
	args = parser.parse_args()

	motif_length = int(args.length.readline())
	sequences = Fasta(args.sequences)

	finder = MotifFinder(sequences, motif_length)

if __name__ == '__main__':
	main()
