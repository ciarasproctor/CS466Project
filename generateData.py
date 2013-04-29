#/usr/bin/python

import os

def main():
	varyICPC()
	varyML()
	varySC()

def varyICPC():
	icpcs = [1.0,1.5,2.0]
	(_,ml,sl,sc) = setDefaults()
	for icpc in icpcs:
		generateDataSets(icpc,ml,sl,sc)

def varyML():
	mls = [6,7,8]
	(icpc,_,sl,sc) = setDefaults()
	for ml in mls:
		generateDataSets(icpc,ml,sl,sc)

def varySC():
	scs = [5.10,20]
	(icpc,ml,sl,_) = setDefaults()
	for sc in scs:
		generateDataSets(icpc,ml,sl,sc)

def setDefaults():
	return (2.0,8,500,10)

def generateDataSets(icpc,ml,sl,sc):
	dirString = '%.1f-%d-%d-%d' % (icpc,ml,sl,sc)
	dirString = os.path.join('data',dirString)
	for i in range(10):
		d = dirString + '_%d' %i
		print(d)
		os.system('./generateBenchmark.py %f %d %d %d %s' %(icpc,ml,sl,sc,d))
		

if __name__ == '__main__':
	main()
