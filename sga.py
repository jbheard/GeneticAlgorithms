'''
TODO:
	- Enforce maximum and minimum values for numerical encoding
	- Don't flip invalid bits for BCD format (e.g. 1001 -> 1101 is invalid, since 1101 is not BCD)
'''

import random
import datetime
import matplotlib.pyplot as plt

from copy import deepcopy
from objfunc import *

# Min and max for floating point numbers
MINRNG, MAXRNG = -5.12, 5.12
NUM_TERMS = 3 # How many 4 bit terms to use for decoding
PRECISION = NUM_TERMS-1 # How many decimal places to add for floating point decoding
BCD = 4 # length of one digit in bits for BCD 

# Get a True value with certain probability, otherwise False
# 0 <= probability <= 1
def flip(probability):
	return random.random() < probability

# Container class for chromosomes
class Genotype:
	def __init__(self, lchrom, fitness=0):
		# Create list of random boolean values
		self.chrom = [ flip(0.5) for _ in range(lchrom) ]
		# Assign fitness value
		self.fitness = fitness
		return

	# For printing out gene values
	def __str__(self):
		ret = '{:.2f} | '.format(self.fitness)
		for x in self.chrom:
			ret += '1' if x else '0'
		return ret

	# evaluate and update fitness
	def eval_fitness(self, obj_func, decode_func):
		self.fitness = obj_func( decode_func(self.chrom) )
		return self.fitness

	def set_chrom(self, newchrom):
		self.chrom = newchrom
		return

###################### Decoding Functions ########################

# Decode a chromosome (in BCD format + sign bit) into signed floating point numbers with given precision
def signed_float_decode(chrom):
	ret = signed_decode(chrom)
	# Get 2 decimal places (divide by 100) for each decoded signed values
	ret = [ x/(10.0**PRECISION) for x in ret ]
	return ret

# Decode a chromosome (in BCD format) into unsigned floating point numbers with given precision
def unsigned_float_decode(chrom):
	ret = unsigned_decode(chrom)
	# Get 2 decimal places (divide by 100) for each decoded signed values
	ret = [ x/(10.0**PRECISION) for x in ret ]
	return ret

# Decode a chromosome (in BCD format + sign bit) into a list of signed integers
def signed_decode(chrom):
	ret = []
	for i in range(0, len(chrom), NUM_TERMS*BCD + 1):
		num = 0
		for j in range(i+1, i+NUM_TERMS*BCD, BCD):
			for k in range(j, j+BCD):
				num += 2**(k-j) if chrom[k] else 0
			num *= 10
		num = num//10 # undo the extra multiply
		if chrom[i]: num = -num # Set negative bit
		ret.append( num ) # Add the num to the list
	return ret

# Decode a chromosome (in BCD format) into a list of unsigned integers
def unsigned_decode(chrom):
	ret = []
	for i in range(0, len(chrom), NUM_TERMS*BCD):
		num = 0
		for j in range(i, i+NUM_TERMS*BCD, BCD):
			for k in range(j, j+BCD):
				num += 2**(k-j) if chrom[k] else 0
			num *= 10
		num = num//10 # undo the extra multiply
		ret.append( num ) # Add the num to the list
	return ret

# Decode a chromosome into 0/1 values
def binary_decode(chrom):
	ret = []
	for x in chrom:
		ret.append(1 if x else 0)
	return ret

# Decode a chromosome into +/- 1 values
def plus_minus_decode(chrom):
	ret = []
	for x in chrom:
		ret.append(1 if x else -1)
	return ret

##################################################################

#################### SGA Operation Functions #####################

# Mutation function, very basic and only changes 1 value in chroma
def mutate(parent, pmutate):
	childgene = parent.chrom[:]

	if flip(pmutate): # Randomly decide whether to mutate
		i = random.randint(0, len(childgene)-1)
		childgene[i] = not childgene[i] # Flip a random gene

	child = Genotype(len(parent.chrom))
	child.set_chrom(childgene)
	return child

# Crossover of two genes
def crossover(parent1, parent2, pcross):
	chrom1, chrom2 = parent1.chrom[:], parent2.chrom[:]
	
	# Cross the two genes
	if flip(pcross):
		i = random.randint(1, len(chrom1)-2)
		chrom1[:i], chrom2[i:] = chrom2[:i], chrom1[i:]
	
	child1 = Genotype(len(chrom1))
	child1.set_chrom(chrom1)
	
	child2 = Genotype(len(chrom2))
	child2.set_chrom(chrom2)

	return child1, child2

# Select() is also technically an operator, but is included as part of the SGA class for ease

##################################################################


class SGA:
	'''
	popsize - Population size
	lchrom - Chromosome length
	pm - mutation probability
	pc - crossover probability
	obj_func - objective function to apply to chromosomes
	minimize - if True, attempt to minimize the solution's obj_func, if False, maximize
	'''
	def __init__(self, popsize, lchrom, pm, pc, obj_func, decode_func, minimize=True):
		self.pc = pc
		self.pm = pm
		self.lchrom = lchrom
		self.popsize = popsize
		self.obj_func = obj_func
		self.decode_func = decode_func
		self.minimize = minimize
		
		# sum and average of fitness for current population
		self.sumfit = 0
		
		# Generate the initial population
		self.pop = [ Genotype(lchrom) for _ in range(popsize) ]
		for gene in self.pop:
			gene.eval_fitness(obj_func, decode_func)
		
		# Select arbitrary best gene and update stats to start
		self.best = self.pop[0]
		self.update_statistics()
		return

	# Select a random gene
	def select(self):
		target = random.uniform(0, self.sumfit)
		fitness = self.sumfit # For minimizing, take complement of each gene fitness
		i = 0
		while i < self.popsize-1 and target > 0:
			if self.minimize:
				fitness -= self.pop[i].fitness
			else:
				fitness = self.pop[i].fitness
			target -= fitness
			i += 1
		return i

	# Evolve the current generation
	def next_generation(self):
		for _ in range(0, self.popsize, 2): # Do this operation self.popsize/2 times
			# Select two genes to cross
			gene1 = self.select()
			gene2 = self.select()
			
			# Cossover our two selected genes
			child1, child2 = crossover(self.pop[gene1], self.pop[gene2], self.pc)
			child1.eval_fitness(self.obj_func, self.decode_func)
			child2.eval_fitness(self.obj_func, self.decode_func)
			self.pop[gene1], self.pop[gene2] = child1, child2
		
		# Mutate all genes with probability pm
		for i in range(self.popsize):
			# Mutate the gene
			self.pop[i] = mutate(self.pop[i], self.pm)
			self.pop[i].eval_fitness(self.obj_func, self.decode_func)
		# Sort the genes in increasing order of fitness
		self.pop.sort(key=lambda g: g.fitness)
		return
	
	# Update information on the current generation
	def update_statistics(self):
		# Find total fitness, best gene for this generation
		self.sumfit = sum(g.fitness for g in self.pop)
		if self.minimize: # Minimize best gene
			gen_best = min(self.pop, key=lambda x: x.fitness) # get best gene of generation
			if gen_best.fitness < self.best.fitness:
				self.best = gen_best
		else: # Maximize best gene
			gen_best = max(self.pop, key=lambda x: x.fitness) # get best gene of generation
			if gen_best.fitness > self.best.fitness:
				self.best = gen_best
		return

if __name__ == '__main__':
	maxgen    =   int(input("Enter maximum # of generations: "))
	popsize   =   int(input("Enter population size ------- > "))
	#lchrom    =   int(input("Enter chromosome length ----- > "))
	pcross    = float(input("Enter crossover probability - > "))
	pmutation = float(input("Enter mutation probability -- > "))
	minimize  =       input("Minimize? (Y/N) ------------- > ")
	minimize  = (minimize.lower() == 'y')

	lchrom = 10 # for objfuncXX, set lchrom = XX
	lchrom = (NUM_TERMS*BCD+1)*2
	#sga = SGA(popsize, lchrom, pmutation, pcross, objfunc10, plus_minus_decode, minimize)
	sga = SGA(popsize, lchrom, pmutation, pcross, himmelblau, signed_float_decode, minimize)
	#sga = SGA(popsize, lchrom, pmutation, pcross, dejong, float_decode, minimize)
	#sga = SGA(popsize, lchrom, pmutation, pcross, rosenbrock, minimize)

	# The main body of the SGA, evolve to next generation and update statistics until max no. of generations have been parsed
	y = []
	x = list(range(maxgen))
	for gen in range(maxgen):
		sga.next_generation()
		sga.update_statistics()
		
		# Add best and average fitness to y vals
		y.append(sga.best.fitness)
	
	print("Best Gene:", str(sga.best) )

	plt.plot(x, y, label='best')
	plt.xlabel('Generation #')
	plt.ylabel('Best, Avg Fitness')
	plt.show()
