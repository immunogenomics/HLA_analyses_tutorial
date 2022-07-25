#!/usr/bin/env python
import sys
missfile = sys.argv[1]
ibdfile = sys.argv[2]
outfile = sys.argv[3]

missing={}
with open(missfile) as f:
	for line in f:
		line=line.rstrip().split()
		if line[0] != "FID":
			missing[(line[0],line[1])]=float(line[5])

REMOVE=open(outfile,"w")

rem_list=[]
with open(ibdfile) as f:
	for line in f:
		line=line.rstrip().split()
		if line[0] != "FID1":
			if missing[(line[0],line[1])] >= missing[(line[2],line[3])]:
				print(line[0]+"\t"+line[1],file=REMOVE)
				rem_list.append(line[0])
			elif missing[(line[0],line[1])] < missing[(line[2],line[3])]:
				print(line[2]+"\t"+line[3],file=REMOVE)
				rem_list.append(line[1])

REMOVE.close()