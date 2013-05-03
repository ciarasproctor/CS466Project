#/usr/bin/python

import os, json

def main():
	directories = os.listdir('./data/')
	sites = {}
	for directory in directories:
		dirString = os.path.join('data',directory)
		siteOverlap = getSiteOverlap(dirString)
		sites[directory] = siteOverlap
	f = open('siteOverlap.json','w')
	json.dump(sites,f,indent=4,sort_keys=True)

def getSiteOverlap(dirName):
	site_file = open(os.path.join(dirName,'sites.txt'))
	predicted_file = open(os.path.join(dirName,'predictedsites_b.txt'),'r')
	sites = readSites(site_file)
	predicted = readSites(predicted_file)
	overlaps = getOverlaps(sites,predicted)
	return overlaps

def readSites(f):
	sites = f.readlines()
	return sites

def getOverlaps(sites,predicted):
	overlaps = []
	for i in range(len(sites)):
		o = int(predicted[i]) - int(sites[i])
		overlaps.append(o)
	return overlaps

if __name__ == '__main__':
	main()
