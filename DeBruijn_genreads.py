class DeBruijnGraph:
   ''' De Bruijn directed multigraph built from a collection of
       strings. User supplies strings and k-mer length k.  Nodes
       are k-1-mers.  An Edge corresponds to the k-mer that joins
       a left k-1-mer to a right k-1-mer. '''

   @staticmethod
   def chop(st, k):
       ''' Chop string into k-mers of given length '''
       for i in range(len(st)-(k-1)):
           yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])

   class Node:
       ''' Node representing a k-1 mer.  Keep track of # of
           incoming/outgoing edges so it's easy to check for
           balanced, semi-balanced. '''

       def __init__(self, km1mer):
           self.km1mer = km1mer
           self.nin = 0
           self.nout = 0

       def isSemiBalanced(self):
           return abs(self.nin - self.nout) == 1

       def isBalanced(self):
           return self.nin == self.nout

       def __hash__(self):
           return hash(self.km1mer)

       def __str__(self):
           return self.km1mer

   def __init__(self, strIter, k, circularize=False):
       ''' Build de Bruijn multigraph given string iterator and k-mer
           length k '''
       self.G = {}     # multimap from nodes to neighbors
       self.nodes = {} # maps k-1-mers to Node objects
       for st in strIter:
           if circularize:
               st += st[:k-1]
           for kmer, km1L, km1R in self.chop(st, k):
               nodeL, nodeR = None, None
               if km1L in self.nodes:
                   nodeL = self.nodes[km1L]
               else:
                   nodeL = self.nodes[km1L] = self.Node(km1L)
               if km1R in self.nodes:
                   nodeR = self.nodes[km1R]
               else:
                   nodeR = self.nodes[km1R] = self.Node(km1R)
               nodeL.nout += 1
               nodeR.nin += 1
               self.G.setdefault(nodeL, []).append(nodeR)
       # Iterate over nodes; tally # balanced, semi-balanced, neither
       self.nsemi, self.nbal, self.nneither = 0, 0, 0
       # Keep track of head and tail nodes in the case of a graph with
       # Eularian walk (not cycle)
       self.head, self.tail = None, None
       for node in iter(self.nodes.values()):
           if node.isBalanced():
               self.nbal += 1
           elif node.isSemiBalanced():
               if node.nin == node.nout + 1:
                   self.tail = node
               if node.nin == node.nout - 1:
                   self.head = node
               self.nsemi += 1
           else:
               self.nneither += 1

   def nnodes(self):
       ''' Return # nodes '''
       return len(self.nodes)

   def nedges(self):
       ''' Return # edges '''
       return len(self.G)



coverage = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
error_rate = [0, 0.01, 0.02, 0.05]
kmer_size = 60
read_length=150

for i in coverage:
   for j in error_rate:
       DeBruijnReads = []
       errorfilename = str(i)+"x_"+ str(read_length) +"nt_"+str(j)+"e.txt"
       error = open(errorfilename, "r")
       for line in error:
           line = line.strip()
           DeBruijnReads.append(line)
       G = DeBruijnGraph(DeBruijnReads, k=kmer_size)
       print(i, j, G.nnodes(), G.nedges(), G.nnodes()+G.nedges())
       error.close()

import random
import math
import numpy as np
import time

def nucleotide2num(char):
   switcher = {
       'A': 0,
       'C': 1,
       'G': 2,
       'T': 3,
   }
   return switcher.get(char)

def num2nucleotide(int):
   switcher = {
       0 : 'A',
       1 : 'C',
       2 : 'G',
       3 : 'T',
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

random.seed(time.time())
coverage = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
read_length = 100
window_range = 48502 - read_length
number_reads = {}


for i in coverage:
   number_reads[i] = int(math.ceil((48502*i)/read_length))

### 0 error rate ###
for i in coverage:
   filename = str(i)+"x_" + str(read_length) + "nt_0e.txt"
   g = open(filename, "w")
   for j in range(number_reads[i]):
       k = random.randint(0,window_range)
       g.write(genome[k:k+read_length])
       g.write("\n")
   g.close()

error_rate = [0.01, 0.02, 0.05]
for i in coverage:
   for j in error_rate:
       errorfilename = str(i)+"x_"+ str(read_length) + "nt_"+str(j)+"e.txt"
       error = open(errorfilename, "w")
       for l in range(number_reads[i]):
           k = random.randint(0,window_range)
           line = genome[k:k+read_length]
           lineerror = mutate(line, j)
           error.write(lineerror)
           error.write("\n")
       error.close()
