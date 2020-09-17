import time
import sys
import csv

bed_file = sys.argv[1]
out_file = sys.argv[2]

nuc = ["A", "T", "G", "C"]

with open(bed_file, "r") as file:
        bed = csv.reader(file, delimiter = "\t")
        for line in bed:
                anc = []
                for i in line[3:9]:
                        if i in nuc:
                                anc.append(i)
                if  len(set(anc)) == 1:
                        out = line[0] + "\t" + line[1] + "\t" + line[2]  + "\t" + anc[0] + "\n"
                        with open(out_file, "a") as target_file:
                                target_file.write(out)
