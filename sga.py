'''
TODO:
	- Instead of just maxgen, add convergence check for optimal gene
'''

import random
import matplotlib.pyplot as plt

from copy import deepcopy
from objfunc import *

# Get a True value with certain probability, otherwise False
# 0 <= probability <= 1
def flip(probability):
	return random.random() < probability

# Flip a random bit for a given digit
def flip_bit_bcd(x):
	x = int(x * (10**PRECISION)) # Convert to integer (keeping precision)
	sign = -1 if x < 0 else 1 # Get sign of number
	
	# Give sign bit a 1/4 chance compared to full digits (BCD)
	k = random.randint(-1, (NUM_DIGITS-1) * 4)

	# Flip sign bit and return
	if k < 0: return -x / (10**PRECISION)
	else: k //= 4
	
	x = abs(x)
	digit = (x // (10**k)) % 10 # Get random digit k
	x -= digit * (10**k)        # Clear digit k

	# Select only values which will result in a valid output
	exclusions = [ j for j in range(0, 9+1) if not (MINRNG <= (x + j*(10**k))/(10**PRECISION) <= MAXRNG) ]

	# Get options for flipping 1 bit based on digit
	if   digit == 0: opts = [1, 2, 4, 8] # 0000 -> 1000, 0100, 0010, 0001
	elif digit == 1: opts = [0, 3, 5, 9] # 1000 -> 0000, 1100, 1010, 1001
	elif digit == 2: opts = [0, 3, 6]    # 0100 -> 1100, 0000, 0110
	elif digit == 3: opts = [1, 2, 7]    # 1100 -> 0100, 1000, 1110
	elif digit == 4: opts = [0, 5, 6]    # 0010 -> 1010, 0110, 0000
	elif digit == 5: opts = [1, 4, 7]    # 1010 -> 0010, 1110, 1000
	elif digit == 6: opts = [2, 4, 7]    # 0110 -> 1110, 0010, 0100
	elif digit == 7: opts = [3, 5, 6]    # 1110 -> 0110, 1010, 1100
	elif digit == 8: opts = [0, 9]       # 0001 -> 1001, 0000
	elif digit == 9: opts = [1, 8]       # 1001 -> 0001, 1000

	# Remove any invalid digits
	opts = [x for x in opts if x not in exclusions]
	if len(opts) == 0: return x
	
	# Set digit k to new digit
	ret = x + random.choice( opts ) * (10**k)
	ret = sign * (ret / (10**PRECISION))
	return ret

def mutate(gene, pmutate, flip_bit):
	updated = False
	if flip(pmutate):
		chrom = gene.chrom
		i = random.randint(0, len(chrom)-1)
		
		chrom[i] = flip_bit(chrom[i])

		# Update chromosome and set updated to true
		gene.set_chrom(chrom)
		updated = True

	return updated

def crossover(parent1, parent2, pcross):
	child1, child2 = deepcopy(parent1), deepcopy(parent2) # Set children

	if flip(pcross) and len(child1.chrom) > 1:
		chrom1, chrom2 = child1.chrom, child2.chrom
		i = random.randint(1, len(chrom1)-1)
		chrom1[:i], chrom2[i:] = chrom2[:i], chrom1[i:]

		# Update chromosomes
		#child1.set_chrom(chrom1)
		#child2.set_chrom(chrom2)

	return child1, child2


# Container class for chromosomes
class Genotype:
	def __init__(self, lchrom, generate):
		# Create list of random boolean values
		self.chrom = [ generate() for _ in range(lchrom) ]
		return

	# For printing out gene values
	def __str__(self):
		ret = '{:.2f} | '.format(self.fitness)
		for x in self.chrom:
			if type(x) is float: # Only print the precision we expect from float
				ret += ('{:.' + str(PRECISION) + 'f} ').format(x)
			else:
				ret += '{} '.format(x)
		return ret

	# evaluate and update fitness
	def eval_fitness(self, obj_func):
		self.fitness = obj_func( self.chrom )
		return self.fitness

	def set_chrom(self, newchrom):
		self.chrom = newchrom
		return

class SGA:
	'''
	popsize - Population size
	lchrom - Chromosome length
	pm - mutation probability
	pc - crossover probability
	obj_func - objective function to apply to chromosomes
	minimize - if True, attempt to minimize the solution's obj_func, if False, maximize
	'''
	def __init__(self, popsize, lchrom, pm, pc, obj_func, generate, flip_bit, minimize=True):
		self.pc = pc
		self.pm = pm
		self.lchrom = lchrom
		self.popsize = popsize
		self.obj_func = obj_func
		self.flip_bit = flip_bit
		self.minimize = minimize
		
		# sum and average of fitness for current population
		self.sumfit = 0
		
		# Generate the initial population
		self.pop = [ Genotype(lchrom, generate) for _ in range(self.popsize) ]
		for gene in self.pop:
			gene.eval_fitness(obj_func)
			
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
			
			child1.eval_fitness(self.obj_func)
			child2.eval_fitness(self.obj_func)

			# Only take children with better fitness than parents
			if child1.fitness > self.pop[gene1].fitness:
				self.pop[gene1] = child1
			if child2.fitness > self.pop[gene2].fitness:
				self.pop[gene2] = child2
		
		# Mutate all genes with probability pm
		for i in range(self.popsize):
			# Mutate the gene
			if mutate(self.pop[i], self.pm, self.flip_bit):
				self.pop[i].eval_fitness(self.obj_func)

		# Sort the genes in increasing order of fitness
		self.pop.sort(key=lambda g: g.fitness)
		return

	# Update information on the current generation
	def update_statistics(self):
		# Find total fitness, best gene for this generation
		self.sumfit = sum(g.fitness for g in self.pop)
		if self.minimize: # Minimize best gene
			gen_best = min(self.pop, key=lambda x: x.fitness) # Get best gene of generation
			if gen_best.fitness < self.best.fitness:
				self.best = deepcopy(gen_best)
		else: # Maximize best gene
			gen_best = max(self.pop, key=lambda x: x.fitness) # Get best gene of generation
			if gen_best.fitness > self.best.fitness:
				self.best = deepcopy(gen_best)

		if self.sumfit > 100000:
			print("sumfit =", self.sumfit)
			for g in self.pop: print(str(g))
			exit(1)
		return

if __name__ == '__main__':
	maxgen    =   int(input("Enter maximum # of generations: "))
	popsize   =   int(input("Enter population size ------- > "))
	#lchrom    =   int(input("Enter chromosome length ----- > "))
	pcross    = float(input("Enter crossover probability - > "))
	pmutation = float(input("Enter mutation probability -- > "))
	minimize  =       input("Minimize? (Y/N) ------------- > ") # if user answers No, maximize instead
	minimize  = (minimize[0].lower() == 'y')

	# GLOBAL VARIABLES
	MINRNG, MAXRNG = -5.00, 5.00 # Min and max for floating point numbers
	NUM_DIGITS = 3 # How many significant digits for each value
	PRECISION = NUM_DIGITS-1 # How many decimal places for floating points (must be <= NUM_DIGITS)

#	lchrom = 27 # for objfuncXX, set lchrom = XX
#	objfunc = objfunc27
#	generate = lambda: random.choice( [-1, 1] ) # lambda for selecting random values
#	sga = SGA(popsize, lchrom, pmutation, pcross, objfunc, generate, lambda x: -x, minimize)

	lchrom = 2
	objfunc = dejong
	generate = lambda: random.uniform( MINRNG, MAXRNG )
	sga = SGA(popsize, lchrom, pmutation, pcross, objfunc, generate, flip_bit_bcd, minimize)

	# The main body of the SGA, evolve to next generation and update statistics until max no. of generations have been parsed
	y = []
	x = list(range(maxgen))
	for gen in range(maxgen):
		sga.next_generation()
		sga.update_statistics()
		
		# TODO: Add min, max, avg fitness to y vals
		y.append(sga.best.fitness)
	
	print("Best Gene:", str(sga.best))
	
	plt.plot(x, y, label='best')
	plt.xlabel('Generation #')
	plt.ylabel('Best, Avg Fitness')
	plt.show()

