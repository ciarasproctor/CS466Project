#/usr/bin/python

import os, json
from pylab import *
import matplotlib.pyplot as plt

def main():
	e_file = open('entropy.json','r')
	entropy = json.load(e_file)
	r_file = open('runningTimes.json','r')
	running = json.load(r_file)
	s_file = open('siteOverlap.json','r')
	sites = json.load(s_file)
	data_sets = generateDataSets(running,entropy,sites)
	print(len(data_sets))
	chartsVaryICPC(data_sets)
	chartsVaryML(data_sets)
	chartsVarySC(data_sets)
	
def chartsVaryICPC(data):
	data_sets = {}
	data_sets[1.0] = getDataSet(1.0,8,10,data)
	data_sets[1.5] = getDataSet(1.5,8,10,data)
	data_sets[2.0] = getDataSet(2.0,8,10,data)
	(x1,r1,s1,e1) = generateArrays(data_sets[1.0],1.0)
	(x15,r15,s15,e15) = generateArrays(data_sets[1.5],1.5)
	(x2,r2,s2,e2) = generateArrays(data_sets[2.0],2.0)
	icpc_arr = [1.0,1.5,2.0]
	arr1 = generateAverages(r1,r15,r2)
	arr2 = generateAverages(s1,s15,s2)
	arr3 = generateAverages(e1,e15,e2)
	plt.xlabel('ICPC')
	plt.ylabel('seconds')
	plt.title('ICPC vs Running Time')
	plt.axis([0,3,0,2])
	plt.plot(x1,r1,'ro',x2,r2,'ro',x15,r15,'ro',icpc_arr,arr1,'g^')
	plt.savefig('runTimeICPC.png')
	plt.clf()
	plt.xlabel('ICPC')
	plt.ylabel('difference in sites')
	plt.title('ICPC vs Site Difference')
	plt.axis([0,3,-2,300])
	plt.plot(x1,s1,'ro',x2,s2,'ro',x15,s15,'ro',icpc_arr,arr2,'g^')
	plt.savefig('sitesICPC.png')
	plt.clf()
	plt.xlabel('ICPC')
	plt.ylabel('entropy')
	plt.title('ICPC vs Entropy')
	plt.axis([0,3,0,1])
	plt.plot(x1,e1,'ro',x2,e2,'ro',x15,e15,'ro',icpc_arr,arr3,'g^')
	plt.savefig('entropyICPC.png')
	plt.clf()
	
def chartsVaryML(data):
	data_sets = {}
	data_sets[6] = getDataSet(2.0,6,10,data)
	data_sets[7] = getDataSet(2.0,7,10,data)
	data_sets[8] = getDataSet(2.0,8,10,data)
	(x1,r1,s1,e1) = generateArrays(data_sets[6],6)
	(x15,r15,s15,e15) = generateArrays(data_sets[7],7)
	(x2,r2,s2,e2) = generateArrays(data_sets[8],8)
	ml_arr = [6,7,8]
	arr1 = generateAverages(r1,r15,r2)
	arr2 = generateAverages(s1,s15,s2)
	arr3 = generateAverages(e1,e15,e2)
	plt.xlabel('ML')
	plt.ylabel('seconds')
	plt.title('ML vs Running Time')
	plt.axis([5,9,0,2])
	plt.plot(x1,r1,'ro',x2,r2,'ro',x15,r15,'ro',ml_arr,arr1,'g^')
	plt.savefig('runTimeML.png')
	plt.clf()
	plt.xlabel('ML')
	plt.ylabel('site differences')
	plt.title('ML vs Site Differences')
	plt.axis([5,9,-2,300])
	plt.plot(x1,s1,'ro',x2,s2,'ro',x15,s15,'ro',ml_arr,arr2,'g^')
	plt.savefig('sitesML.png')
	plt.clf()
	plt.xlabel('ML')
	plt.ylabel('entropy')
	plt.title('ML vs Entropy')
	plt.axis([5,9,0,1])
	plt.plot(x1,e1,'ro',x2,e2,'ro',x15,e15,'ro',ml_arr,arr3,'g^')
	plt.savefig('entropyML.png')
	plt.clf()

def chartsVarySC(data):
	data_sets = {}
	data_sets[5] = getDataSet(2.0,8,5,data)
	data_sets[10] = getDataSet(2.0,8,10,data)
	data_sets[20] = getDataSet(2.0,8,20,data)
	(x1,r1,s1,e1) = generateArrays(data_sets[5],5)
	(x15,r15,s15,e15) = generateArrays(data_sets[10],10)
	(x2,r2,s2,e2) = generateArrays(data_sets[20],20)
	sc_arr = [5,10,20]
	arr1 = generateAverages(r1,r15,r2)
	arr2 = generateAverages(s1,s15,s2)
	arr3 = generateAverages(e1,e15,e2)
	plt.xlim(-1,25)
	plt.xlabel('SC')
	plt.ylabel('seconds')
	plt.title('SC vs Running Time')
	plt.plot(x1,r1,'ro',x2,r2,'ro',x15,r15,'ro',sc_arr,arr1,'g^')
	plt.savefig('runTimeSC.png')
	plt.clf()
	plt.xlim(-1,25)
	plt.xlabel('SC')
	plt.ylabel('site differences')
	plt.title('SC vs Site Differences')
	plt.xlim(-1,25)
	plt.ylim(-1,250)
	plt.plot(x1,s1,'ro',x2,s2,'ro',x15,s15,'ro',sc_arr,arr2,'g^')
	plt.savefig('sitesSC.png')
	plt.clf()
	plt.xlabel('SC')
	plt.ylabel('entropy')
	plt.title('SC vs Entropy')
	plt.xlim(-1,25)
	plt.ylim(0,1)
	plt.plot(x1,e1,'ro',x2,e2,'ro',x15,e15,'ro',sc_arr,arr3,'g^')
	plt.savefig('entropySC.png')
	plt.clf()

def generateArrays(data,x):
	x_arr = []
	r_arr = []
	s_arr = []
	e_arr = []
	for d in data:
		x_arr.append(x)
		r_arr.append(d.run_time)
		s_arr.append(d.avgSite())
		e_arr.append(d.avgEntropy())
	return(x_arr,r_arr,s_arr,e_arr)

def generateAverages(arr1,arr2,arr3):
	r = avgList(arr1)
	s = avgList(arr2)
	e = avgList(arr3)
	return [r,s,e]

def avgList(l):
	s = 0
	for item in l:
		s += item
	return float(s)/len(l)


def generateDataSets(running,entropy,sites):
	data_sets = []
	for data in entropy.keys():
		spl = data.split('-')
		i = float(spl[0])
		m = int(spl[1])
		s = int(spl[3][:spl[3].index('_')])
		ds = DataSet(i,m,s,running[data],entropy[data],sites[data])
		data_sets.append(ds)
	return data_sets

def getDataSet(i,m,s,data):
	data_set = []
	for d in data:
		if d.icpc == i and d.ml == m and d.sc == s:
			data_set.append(d)
	return data_set	

class DataSet:
	def __init__(self,icpc,ml,sc,rt,e,s):
		self.icpc = icpc
		self.ml = ml
		self.sc = sc
		self.run_time = rt
		self.entropy = e
		self.site = s
	
	def avgEntropy(self):
		s = 0
		print(self.entropy)
		for e in self.entropy:
			s+= e
		return float(s)/len(self.entropy)

	def avgSite(self):
		s = 0
		for site in self.site:
			s+= abs(site)
		return float(s)/len(self.site)

if __name__ == '__main__':
	main()
