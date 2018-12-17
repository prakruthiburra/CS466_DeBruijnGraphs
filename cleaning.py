f = open("Raw_genome.txt")
sequence = ""
for line in f:
    line = line.split()
    for i in range(1,len(line)):
        sequence = sequence+line[i]
f.close()

f = open("Cleaned_genome.txt", "w")
f.write(sequence)

print(len(sequence))
