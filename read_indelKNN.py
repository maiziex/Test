#!/usr/bin/env python
import pdb 
import matplotlib.pyplot as plt
import pylab
import numpy as np
import os
import sys


def load_indel_knn():
    indel_knn = {}
     
    with open('Software/population_INDEL/phase3_result_cutblock.final', 'r') as inFile:
        file = inFile.read()
             
    inFile.close()
     
    for line in file.strip().split('\n'):
        try:
        	[locus,indel_id,indel_len,af,SI,locus1,si1,locus2,si2,locus3,si3,locus4,si4,empty] = line.split('\t')
        	indel_knn[locus] = [indel_id,indel_len,af,SI,locus1,si1,locus2,si2,locus3,si3,locus4,si4]
	except ValueError:
		[locus,indel_id,indel_len,af,SI,locus1,si1,locus2,si2,locus3,si3,locus4,si4] = line.split('\t')
		indel_knn[locus] = [indel_id,indel_len,af,SI,locus1,si1,locus2,si2,locus3,si3,locus4,si4]
 
    return indel_knn
     
 
def main():
    pdb.set_trace()
    count_noused = 0
    count_used = 0
    indel_knn = load_indel_knn()
    data = {"indel_len":[],"x":[], "y":[], "z":[],"label":[]}
    for label, coord in indel_knn.items(): 
	if  int(coord[0]) == 1:                
		if float(coord[2])<0.05:
		#if (abs(int(coord[4])- int(label))+abs(int(coord[6])- int(label))+abs(int(coord[8])- int(label))+abs(int(coord[10])- int(label)))/float(4)>800000:	
			count_used = count_used+1
        		data["indel_len"].append(float(coord[1]))
 			data["x"].append(float(coord[2])) 		
        		data["y"].append((abs(int(coord[4])- int(label))+abs(int(coord[6])- int(label))+abs(int(coord[8])- int(label))+abs(int(coord[10])- int(label)))/float(4))
			#data["z"].append((float(coord[5])+float(coord[7])+float(coord[9])+float(coord[11]))/float(4))
			data["z"].append(float(coord[3]))
    			data["label"].append(abs(int(coord[0])))
                        print label, coord[0],coord[1],coord[2],coord[3],coord[4],coord[6],coord[8],coord[10],coord[5],coord[7],coord[9],coord[11]
		else:
			count_noused = count_noused+1
		
    x = data["x"]
    y = data["z"] 
    indel_len = data["indel_len"]
    # plot the data itself
    #pylab.plot(x,y,'.')
    #numBins=50
    n, bins, patches = plt.hist(indel_len, 50, normed=1, facecolor='green', alpha=0.75)
    plt.ylabel('strength')
    plt.xlabel('af')
    print count_noused, count_used
    plt.show()

 
    return
 
if __name__ == '__main__':
    main()
