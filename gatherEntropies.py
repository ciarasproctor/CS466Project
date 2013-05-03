#/usr/bin/python

import os, math, json

def main():
	directories = os.listdir('./data/')
	data = {}
	for directory in directories:
		dirString = os.path.join('data',directory)
		entropy = getEntropy(dirString)
		data[directory] = entropy
	f = open('entropy.json','w')
	json.dump(data,f,indent=4,sort_keys=True)

def getEntropy(dirName):
	#open files
	motif_file = open(os.path.join(dirName,'motif.txt'),'r')
	predicted_file = open(os.path.join(dirName,'predictedmotif_b.txt'),'r')
	#read matrices
	motif = readMotif(motif_file)
	predicted = readMotif(predicted_file)
	#calculate entropy
	entropy = getPositionalEntropy(motif,predicted)
	#return
	return entropy

def readMotif(f):
	motif = []
	for line in f.readlines()[1:]:
		if line[0] != '<':
			sp = line.split()
			f_sp = []
			for s in sp:
				f_sp.append(float(s))
			motif.append(f_sp)
	return motif

def getPositionalEntropy(motif,predicted):
	pe = []
	print(motif)
	print(motif[0][0], motif[0][1], motif[0][2], motif[0][3])
	sc = motif[0][0] + motif[0][1] + motif[0][2] + motif[0][3] 
	for p in range(len(motif)):
		e = 0.0
		for b in range(4):
			if motif[p][b] != 0 and predicted[p][b] != 0:
				e += (motif[p][b] / sc)* math.log(((motif[p][b]/sc)/(predicted[p][b]/sc)),2)
		e = e/len(motif)	
		pe.append(e)
	return pe

if __name__ == '__main__':
	main()
