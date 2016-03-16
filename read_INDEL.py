#!/usr/bin/env python
import pdb 
import matplotlib.pyplot as plt
import pylab
import numpy as np
import os
import sys
from operator import itemgetter
#pdb.set_trace()

  

def load_indel_all():
    indel_all = {}
     
    with open('Software/population_INDEL/saved_indel.txt', 'r') as inFile:
        file = inFile.read()
             
    inFile.close()
     
    for line in file.strip().split('\n'):
        try:
        	[locusandindel,count] = line.split('\t')
        	locus = locusandindel[0:8]
	        indel = locusandindel[8:]
		if locus in indel_all.keys():
                	indel_all[locus].append([indel,count])
		else:
			indel_all[locus] = [[indel,count]]
	except ValueError:
		print "wrong line"
		
 
    return indel_all
        
 
def main():
    pdb.set_trace()
    indel_all = load_indel_all()
    data = {"x":[], "y":[]}
    all_lenofalt = []
    f = open('indel_vcf_raw','w')
    for label, coord in indel_all.items():
	f.write(label)
	f.write("\t")	
    	for ii in range(0,len(coord)):			
        	data["x"].append(int(coord[ii][1]))
		data["y"].append(coord[ii][0])
	
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
	data = {"x":[], "y":[]}	
	

    f.close()	
    print all_lenofalt
  
 
    return
 
if __name__ == '__main__':
    main()

