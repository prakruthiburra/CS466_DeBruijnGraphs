import timeit

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

kmer_size=30
l=range(2500, 48502, 2500)

filename = "Cleaned_genome.txt"
f = open(filename, "r")
genome = f.read().replace('\n', '')
f.close()

for temp in l:
   start = timeit.default_timer()
   DeBruijnGraph(genome[:temp], k=kmer_size)
   stop = timeit.default_timer()
   print('Length:', temp, ',', 'Time: ', stop - start)
