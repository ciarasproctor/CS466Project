#!/usr/bin/python
import argparse, random, numpy, math, os

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
	parser.add_argument('directory', type=str, help='Write output files to this directory', default='.')
	args = parser.parse_args()

	createSubdirectory(args.directory)
	sequences = generateSequences(args.SC, args.SL)
	motif = generateMotif(args.ML, args.ICPC, args.SC)
	bindingSites = generateBindingSites(motif,args.SC)
	sites = plantBindingSites(sequences,bindingSites, args.SL, args.ML)
	writeSequences(sequences,args.directory)
	writeSites(sites,args.directory)
	writeMotif(motif,args.ML,args.directory)
	writeMotifLength(args.ML,args.directory)

def createSubdirectory(dirName):
	if not os.path.exists(dirName):
		os.makedirs(dirName)

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

def generateMotif(ml, icpc, sc):
	if icpc == 2:
		return generateMotifIC2(ml,sc)
	target = ml*icpc*-1.0
	pwm = numpy.zeros(shape=(ml, 4))
	for row in range(ml):
		generateMotifColumns(pwm, row, sc)
	#fix information content
	ic = getInformationContent(pwm,sc)
	while ic != target:
		pwm = adjustInformationContent(pwm,sc,ic,target)
		ic = getInformationContent(pwm,sc)
		print('ic: %.20f' % ic)
	return pwm

def generateMotifIC2(ml,sc):
	pwm = numpy.zeros(shape=(ml,4))
	for row in range(ml):
		generateMotifIC2Columns(pwm,row,sc)
	return pwm

def generateMotifIC2Columns(pwm,row,sc):
	colNum = random.randint(0,3)
	pwm[row,colNum] = sc

def changeInIC(sc,c1,c2,x):
	if x == 0:
		return 0
	newIC = getInformationContentBase(c1+x,sc) + getInformationContentBase(c2-x,sc)
	oldIC = getInformationContentBase(c1,sc) + getInformationContentBase(c2,sc)
	return abs(newIC) - abs(oldIC)

def adjustInformationContent(pwm,sc,ic,target):
	bound = target-ic
	if abs(ic) < abs(target):
		return raiseInformationContent(pwm,sc,bound)
	else:
		return lowerInformationContent(pwm,sc,bound)

def getCellsChange(pwm):
	rowNum = random.randint(0,len(pwm)-1)
	cols = [0,1,2,3]
	colNum1 = random.choice(cols)
	cols.remove(colNum1)
	colNum2 = random.choice(cols)
	changeValue = random.randint(0,pwm[rowNum,colNum2])
	return (rowNum,colNum1,colNum2,changeValue)

def raiseInformationContent(pwm,sc,bound):
	(rowNum,colNum1,colNum2,change) = getCellsChange(pwm)
	icChange = changeInIC(sc,pwm[rowNum,colNum1],pwm[rowNum,colNum2],change)
	while icChange <= 0 and icChange > bound:
		(rowNum,colNum1,colNum2,change) = getCellsChange(pwm)
		icChange = changeInIC(sc,pwm[rowNum,colNum1],pwm[rowNum,colNum2],change)
	pwm[rowNum,colNum1] = pwm[rowNum,colNum1] + change
	pwm[rowNum,colNum2] = pwm[rowNum,colNum2] - change
	print('raised by: %.20f' % icChange)
	return pwm
	
def lowerInformationContent(pwm,sc,bound):
	(rowNum,colNum1,colNum2,change) = getCellsChange(pwm)
	icChange = changeInIC(sc,pwm[rowNum,colNum1],pwm[rowNum,colNum2],change)
	while icChange >= 0 and icChange < bound:
		(rowNum,colNum1,colNum2,change) = getCellsChange(pwm)
		icChange = changeInIC(sc,pwm[rowNum,colNum1],pwm[rowNum,colNum2],change)
	pwm[rowNum,colNum1] = pwm[rowNum,colNum1] + change
	pwm[rowNum,colNum2] = pwm[rowNum,colNum2] - change
	print('lowered by: %.20f' % icChange)
	return pwm	

def generateMotifColumns(pwm, row, sc):
	samples = sc
	nucleotides = [0,1,2,3]
	for _ in range(3):
		col = random.choice(nucleotides)
		nucleotides.remove(col)
		frequency = random.randint(0,samples)
		pwm[row,col] = frequency
		samples = samples - frequency
	pwm[row,random.choice(nucleotides)] = samples

def getInformationContent(pwm, sc):
	ic = 0
	for (x,y), value in numpy.ndenumerate(pwm):
		if value != 0:
			ic = ic + getInformationContentBase(value,sc)
	return -1*ic

def getInformationContentBase(value,sc):
	if value == 0:
		return 0
	return value/sc * math.log((value/sc)/.25,2)

def generateBindingSites(motif,sc):
	bindingSites = []
	for _ in range(sc):
		bindingSites.append(generateBindingSite(motif,sc))
	return bindingSites

def generateBindingSite(motif,sc):
	bindingSite = ''
	bases = ['A','C','G','T']
	for row in motif:
		count = sc
		baseSelector = random.randint(1,count)
		for i in range(len(row)):
			column = row[i]
			baseSelector -= column
			if baseSelector <= 0:
				bindingSite += bases[i]
				break
	return bindingSite

def plantBindingSites(sequences,bindings,sl,ml):
	plantStartMax = sl-ml
	sites = []
	for i in range(len(sequences)):
		startSite = random.randint(0,plantStartMax)
		sq = sequences[i]
		sequences[i] = sq[:startSite] + bindings[i] + sq[startSite+ml:]
		sites.append(startSite)
	return sites

def writeSequences(sequences,prefix):
	f = open(os.path.join(prefix,'sequences.fa'), 'w')
	lines = []
	for i in range(len(sequences)):
		f.write('>%d \n' % i)
		f.write(sequences[i]+'\n')
	f.writelines(lines)
	f.close()

def writeSites(sites,prefix):
	f = open(os.path.join(prefix,'sites.txt'),'w')
	for site in sites:
		f.write(str(site)+'\n')
	f.close()

def writeMotif(pwm,ml,prefix):
	f = open(os.path.join(prefix,'motif.txt'),'w')
	f.write('>motif %d \n' %ml)
	for row in pwm:
		rowText = ""
		for cell in row:
			rowText += str(cell) + " "
		rowText=rowText.strip()
		f.write(rowText+'\n')
	f.write('<\n')
	f.close()

def writeMotifLength(ml,prefix):
	f = open(os.path.join(prefix,'motiflength.txt'),'w')
	f.write('%s\n' % str(ml))
	f.close()

if __name__ == '__main__':
	main()
