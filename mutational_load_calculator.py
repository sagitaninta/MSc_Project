#### Import modules ####
import sys
import csv
import time

start = time.time()

##### Functions #####




#### Code ####

phylop_file = sys.argv[1]
out_file = sys.argv[2]

all_scores = 0
homozygous_positions = 0

with open(phylop_file, "r") as file:
	phylop_scores = csv.reader(file, delimiter = "\t")
	for line in phylop_scores:
		if line[7] in line[4]:
			score = (1 - float(line[3]))*float(line[6])
		elif line[7] not in line[4]:
			score = float(line[3])*float(line[6])
		all_scores += score
		homozygous_positions += 1

print(all_scores)
print(homozygous_positions)
mutational_load = all_scores/homozygous_positions
print(mutational_load)

#with open(out_file, "a") as file:
#	file.write("Raw PhyloPi score: " + "\t" + str(all_scores) + "\n")
#	file.write("Total homozygous positions: " + "\t" + str(homozygous_positions) + "\n")
#	file.write("Mutational load across homozygous positions:" + "\t" + str(mutational_load))


with open(out_file, "a") as file:
	file.write("PhyloPi score" + "\t" + "Homozygous positions" + "\t" + "Load" + "\n")
	file.write(str(all_scores) + "\t" + str(homozygous_positions) + "\t" + str(mutational_load) + "\n")
