#!/usr/bin/env python
import gzip
import sys

INPUTVCF = sys.argv[1]
CONVERTER = sys.argv[2]
OUTINFOMATRIX =sys.argv[3]
OUT = open(OUTINFOMATRIX, "w")
out = "\t".join(["SNP","CHRPOS","REF","ALT","R2","AF"])
print(out, file = OUT)

variant_dic = {}
with open(CONVERTER, "r") as f:
	for line in f:
		line = line.rstrip().split()
		POS = line[0]
		SNP = line[1]
		A1 = line[2]
		A2 = line[3]
		variant_dic.setdefault(":".join(["6",POS]), []).append([SNP,A1,A2])

with gzip.open(INPUTVCF, "rt", "utf_8") as f:
	for line in f:
		line = line.rstrip()
		if line.find("#") > -1:
			print(line)
		else:
			line = line.split("\t")
			variant =line[2]
			info = line[7].split(";")
			ref = line[3]
			alt = line[4]
			converted = "NA"
			for cand in variant_dic[variant]:
				if (ref == cand[1] and alt == cand[2]) or (ref == cand[2] and alt == cand[1]):
					converted = cand[0]
			out = "\t".join([line[0],line[1],converted]) + "\t" + "\t".join(line[3:])
			print(out)
			info = line[7].split(";")
			for i in range(len(info)):
				if info[i].find("R2") == 0:
					R2 = info[i].split("=")[1]
				elif info[i].find("AF") == 0:
					AF = info[i].split("=")[1]
			out = "\t".join([converted,variant,ref,alt,R2,AF])
			print(out, file = OUT)
OUT.close()

