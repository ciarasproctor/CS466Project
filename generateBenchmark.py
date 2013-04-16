import sys, random

def main():
	if len(sys.argv) < 5:
		print("usage: generateBenchmark.py <ICPC> <ML> <SL> <SC>")
	icpc = float(sys.argv[1])
	ml = int(sys.argv[2])
	sl = int(sys.argv[3])
	sc = int(sys.argv[4])
	sequences = generateSequences(sc, sl)
	
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
