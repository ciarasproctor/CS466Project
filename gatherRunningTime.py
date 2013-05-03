#/usr/bin/python

import os, json

def main():
	directories = os.listdir('./data/')
	times = {}
	for directory in directories:
		dirString = os.path.join('data',directory)
		rt = getRunningTime(dirString)
		times[directory] = rt
	f = open('runningTimes.json','w')
	json.dump(times,f,indent=4,sort_keys=True)

def getRunningTime(dirName):
	time_file = open(os.path.join(dirName,'running_time_b.txt'),'r')
	time = float(time_file.readline())
	return time

if __name__ == '__main__':
	main()
