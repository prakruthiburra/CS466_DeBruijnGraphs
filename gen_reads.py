import random
import math
import numpy as np

def nucleotide2num(char):
    switcher = {
        'a': 0,
        'c': 1,
        'g': 2,
	't': 3,
    }
    return switcher.get(char)

def num2nucleotide(int):
    switcher = {
        0 : 'a',
        1 : 'c',
        2 : 'g',
        3 : 't',
    }
    return switcher.get(int)

def mutate(orline, er):
    newline = list(orline)
    for i in range(0, read_length):
        k = random.random()
        if(k<=er):
            originalchar = [nucleotide2num(orline[i])]
            num = [0, 1, 2, 3]
            choicechar = np.setdiff1d(num,originalchar, assume_unique=True)
            c = random.choice(choicechar)
            newline[i] = num2nucleotide(c)
    newline = ''.join(newline)
    return newline

f = open("Cleaned_genome.txt", "r")
genome = f.read().replace('\n', '')
f.close()

random.seed(1)
coverage = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
read_length = 100
window_range = 48502 - read_length
number_reads = {}


for i in coverage:
    number_reads[i] = int(math.ceil((48502*i)/read_length))

### 0 error rate ###
for i in coverage:
    filename = str(i)+"x_100nt_0e.txt"
    g = open(filename, "w")
    for j in range(number_reads[i]):
        k = random.randint(0,window_range)
        g.write(genome[k:k+read_length])
        g.write("\n")
    g.close()

error_rate = [0.01, 0.02, 0.05]
for i in coverage:
    for j in error_rate:
        count = 0
        errorfilename = str(i)+"x_100nt_"+str(j)+"e.txt"
        error = open(errorfilename, "w")
        for l in range(number_reads[i]):
            count += 1
            k = random.randint(0,window_range)
            line = genome[k:k+read_length]
            print(line)
            print(i, j, count)
            lineerror = mutate(line, j)
            error.write(lineerror)
            error.write("\n")
        error.close()
