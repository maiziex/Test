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
    data = {"x":[], "y":[],"z":[],"yz":[]}
    dict1={}
    dict2={}
    dict3={}
    all_lenofalt = []
    f = open('indel_vcf_addcriterion_raw','w')
    for label, coord in indel_all.items():
	f.write(label)
	f.write("\t")	
    	for ii in range(0,len(coord)):			
        	data["x"].append(int(coord[ii][1]))    #count
		data["y"].append(coord[ii][0])   #indel
		data["z"].append(int(coord[ii][2]))  # samplenum
		data["yz"].append([coord[ii][0],int(coord[ii][2])])  # indel and samplenum
	
	for jj in range(len(data["x"])):
		if data["y"][jj] in dict1.keys():
			dict1[data["y"][jj]]+= data["x"][jj]
		else:
			dict1[data["y"][jj]] = data["x"][jj]

	for kk in range(len(data["x"])):
		if data["yz"][kk] in dict2.keys():
			dict2[str(data["yz"][kk])].append(data["x"][kk])
		else:
			dict2[str(data["yz"][kk])] = [data["x"][kk]]

	for nn in range(len(data["x"])):
		if data["y"][nn] in dict3.keys():
			dict3[str(data["y"][nn])].append([data["z"][nn]])
		else:
			dict3[str(data["y"][nn])] = [[data["z"][nn]]]


	#if len(dict1)>2:
	#	print "here"

        d_view = [ (v,k) for k,v in dict1.iteritems() ]
	d_view.sort(reverse=True) # natively sort tuples by first element
	#for v,k in d_view:
    		#print("There are " + str(v) + " " + k)
	
       
       
        alt_count = 1
        countp = 1
	for v,k in d_view:
		if countp < 3:
			#print label, data["y"][ll], data["x"][ll]		
			f.write(k)
			f.write("\t")
			f.write(str(v))
			f.write("\t")
		else:
			#print "here"	
			if len(dict3[k]) > 1:    # two individual has the variant
				for mm in range(len(dict3[k])):			             
					count = dict2[str([k]+dict3[k][mm])]
					if count[0] > 1:      # at least 2 reads for each individual
						count = [0]
						continue
					else:
						break
						                 
					f.write(k)
					f.write("\t")
					f.write(str(v))
					f.write("\t")
					alt_count = alt_count + 1
					


		countp =  countp +1
		
	all_lenofalt.append(alt_count)		
 	f.write(os.linesep)

	
	data = {"x":[], "y":[],"z":[],"yz":[]}	
	dict1={}
        dict2={}
	dict3={}
	d_view =[]

    f.close()	
    print all_lenofalt
  
 
    return
 
if __name__ == '__main__':
    main()

