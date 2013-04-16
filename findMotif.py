#!/usr/bin/python

import argparse, random

def main():
	parser = argparse.ArgumentParser(description='Motif finding using the Gibbs Motif Sampling Algorithm', formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='''Generates the following files (optionally prefixed by prefix):
  * predictedmotif.txt
  * predictedsites.txt''')
	parser.add_argument('sequences', type=file, help='The file containing the sequences (FASTA format)')
	parser.add_argument('length', type=file, help='Path to the file containing the motif length')
	parser.add_argument('--prefix', type=str, help='Prefix output files with this string', default='')
	args = parser.parse_args()

if __name__ == '__main__':
	main()
