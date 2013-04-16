#!/usr/bin/python
import argparse, random

def main():
	parser = argparse.ArgumentParser(description='Motif benchmark generator', formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='''Generates the following files (optionally prefixed by prefix):
  * sequences.fa
  * sites.txt
  * motif.txt
  * motiflength.txt''')
	parser.add_argument('ICPC', type=int, help='The information content per column')
	parser.add_argument('ML', type=int, help='The motif length')
	parser.add_argument('SL', type=int, help='The length of the generated sequences')
	parser.add_argument('SC', type=int, help='The number of sequences to generate')
	parser.add_argument('--prefix', type=str, help='Prefix output files with this string', default='')
	args = parser.parse_args()

	sequences = generateSequences(args.SC, args.SL)

def generateSequence(sl):
	seq = ""
	for _ in range(sl):
		seq += random.choice("ACTG")
	return seq

def generateSequences(sc, sl):
	sequences = []
	for _ in range(sc):
		sequences.append(generateSequence(sl))
	return sequences


if __name__ == '__main__':
	main()
