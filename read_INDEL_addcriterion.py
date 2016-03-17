#!/usr/bin/env python
import pdb 
import matplotlib.pyplot as plt
import pylab
import numpy as np
import os
import sys
from operator import itemgetter
#pdb.set_trace()

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False  

def load_indel_all():
    indel_all = {}
     
    with open('Software/population_INDEL/saved_indel_addcriterion.txt', 'r') as inFile:
        file = inFile.read()
             
    inFile.close()
    
    for line in file.strip().split('\n'):
	locus  = ""
	indel  = ""
	samplenum = ""
        try:
        	[locusandindel,count] = line.split('\t')
		for ii in range(len(locusandindel)):
			if is_number(locusandindel[ii]): 
        			locus = locus + locusandindel[ii]
			else:
				indel = indel + locusandindel[ii]				
				break
		for jj in range(ii+1,len(locusandindel)):
			if is_number(locusandindel[jj]) and locusandindel[jj-1] != "+" and locusandindel[jj-1] != "-": 
				samplenum = samplenum + locusandindel[jj]
			else:
				indel = indel + locusandindel[jj]
			
		if locus in indel_all.keys():
                	indel_all[locus].append([indel,count,samplenum])
		else:
			indel_all[locus] = [[indel,count,samplenum]]
	except ValueError:
		print "wrong line"
		
 
    return indel_all
        
 
def main():
    pdb.set_trace()
    indel_all = load_indel_all()
    data = {"x":[], "y":[],"z":[]}
    all_lenofalt = []
    f = open('indel_vcf_addcriterion_raw','w')
    for label, coord in indel_all.items():
	f.write(label)
	f.write("\t")	
    	for ii in range(0,len(coord)):			
        	data["x"].append(int(coord[ii][1]))    #count
		data["y"].append(coord[ii][0])   #indel
		data["z"].append(int(coord[ii][2]))  # samplenum
	
	indices, L_sorted = zip(*sorted(enumerate(data["x"]), key=itemgetter(1)))
        lenofalt = len(indices)-1
        all_lenofalt.append(lenofalt)
        for jj in reversed(indices):
		#print label, data["y"][jj], data["x"][jj]		
		f.write(data["y"][jj])
		f.write("\t")
		f.write(str(data["x"][jj]))
		f.write("\t")		
 	f.write(os.linesep)	
	data = {"x":[], "y":[],"z":[]}	
	

    f.close()	
    print all_lenofalt
  
 
    return
 
if __name__ == '__main__':
    main()

