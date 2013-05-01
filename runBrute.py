#/usr/bin/python

import os

def main():
	directories = os.listdir('./data/')
	for directory in directories:
		dirString = os.path.join('data',directory)
		print('running on %s' %dirString)
		os.system('./findMotifBrute.py %s' % dirString)

if __name__ == '__main__':
	main()
