import subprocess, os, string
import multiprocessing
import csv
import math
import getopt
import sys
import os.path
import argparse
import re
import gzip
import random
from functools import partial
from collections import OrderedDict


parser = argparse.ArgumentParser(description="Mutate given sequence")
parser.add_argument("-s", "--sequence", default="ATCGAGCTACCGATTAGGCCGATA",	type=str, help="-s [--sequence] Sequence to mutate")
parser.add_argument("-a", "--allSingleBP", default="T", type=bool, help="-a [--allSingleBP] If TRUE output all single base mutations")
parser.add_argument("-m", "--makeWindow", default="T", type=bool, help="-m [--makeWindow] If TRUE will output windowed mutations, and apply options -w and -s below")
parser.add_argument("-w", "--window", default="8", type=int, help="window size to mutate")
parser.add_argument("-st", "--step", default="4", type=int, help="step size for mutations, recomend ~1/2 of window size")
args = parser.parse_args()

seq = args.sequence
allSingleBP = args.allSingleBP
makeWindow = args.makeWindow
window = args.window
step = args.step

allNucs = ["A","C","G","T"]
pure = ["A","G"]
pyrim = ["T","C"]


def randomSingleMutation(nuc):
	allNucs = ["A","C","G","T"]
	newNucs = allNucs
	newNucs.remove(nuc)
	random.shuffle(newNucs)
	randomNuc = newNucs[0]
	return random

def singleOrderedMutation(nuc,n):
	allNucs = ["A","C","G","T"]
	newNucs = allNucs
	newNucs.remove(nuc)
	return newNucs[n]

def randomTransition(nuc):
	if nuc in pure:
		newNucs = pyrim
		random.shuffle(newNucs)
		randomNuc = newNucs[0]
		return randomNuc

	if nuc in pyrim:
		newNucs = pure
		random.shuffle(newNucs)
		randomNuc = newNucs[0]
		return randomNuc

def makeAllSingleMutations(input):
	allSeqsList = [] # make empty list to file
	allSeqsList.append(input) # add reference seq
	counter = 0
	lengthSeq = len(input)
	for bp in input:
		start = counter + 1
		previousStop = counter
		if previousStop < 0:
			previousStop = 0
		newStartSeq = input[0:previousStop]
		newEndSeq = input[start:lengthSeq]
		for i in range(0,3):
			mutBP = singleOrderedMutation(bp,i)
			newSeq = str(newStartSeq) + str(mutBP) + str(newEndSeq)
			allSeqsList.append(newSeq)
		counter += 1
	return allSeqsList


def makeWindowedTransverstions(input,windowN,stepN):
	allTransList = []
	counter = 0
	lengthSeq = len(input)
	for startLoc in range(0, lengthSeq, stepN):
		stopLoc = startLoc + windowN
		chunk =  input[startLoc:stopLoc]
		previousStop = startLoc
		newStartSeq = input[0:previousStop]
		newEndSeq = input[stopLoc:lengthSeq]
		chunkList = []
		for i in range(0,8):
			newChunk = ''
			for bp in chunk:
				newBP = randomTransition(bp)
				newChunk = newChunk + newBP
			while newChunk in chunkList:
				newBP = randomTransition(bp)
				newChunk = newChunk + newBP
			if newChunk not in chunkList:
				chunkList.append(newChunk)
			newSeq = str(newStartSeq) + str(newChunk) + str(newEndSeq)
			allTransList.append(newSeq)

	return allTransList


seqList = makeAllSingleMutations(seq)

secondSeqList = makeWindowedTransverstions(seq,window,step)

counter=0
for j in secondSeqList:
	#counter+=1
	#print ">Single_mut_"+str(counter)
	print(j)





