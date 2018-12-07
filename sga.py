import random
from copy import deepcopy

TOURN_MAX = 3

# Get a True value with certain probability, otherwise False
# 0 <= probability <= 1
def flip(probability):
	return random.random() < probability

def mutate(gene, pmutate, flip_bit):
	updated = False
	if flip(pmutate):
		flip_bit(gene.chrom)
		# Update chromosome and set updated to true
		updated = True

	return updated

def crossover(parent1, parent2, pcross):
	child1, child2 = deepcopy(parent1), deepcopy(parent2) # Set children

	if flip(pcross) and len(child1.chrom) > 1:
		chrom1, chrom2 = child1.chrom, child2.chrom
		i = random.randint(1, len(chrom1)-1)
		chrom1[:i], chrom2[:i] = chrom2[:i], chrom1[:i]

	return child1, child2

# Container class for chromosomes
class Genotype:
	def __init__(self, lchrom, generate):
		# Create list of random boolean values
		self.chrom = generate(lchrom)
		return

	# For printing out gene values
	def __str__(self):
		ret = '{:.2f} | '.format(self.fitness)
		for x in self.chrom[:-1]:
			if type(x) is float: # Format floats to drop trailing 0s
				ret += '{:.3f},'.format(x)
			else:
				ret += '{:2},'.format(x)
		if type(x) is float: # Format floats to drop trailing 0s
			ret += '{:.3f}'.format(x)
		else:
			ret += '{:2}'.format(x)
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
	def __init__(self, popsize, lchrom, pm, pc, obj_func, generate, flip_bit, minimize=True, cross=crossover, mutate=mutate):
		self.pc = pc
		self.pm = pm
		self.lchrom = lchrom
		self.popsize = popsize
		self.flip_bit = flip_bit
		self.minimize = minimize
		
		# Functions
		self.obj_func  = obj_func
		self.crossover = cross
		self.mutate    = mutate
		
		# sum and average of fitness for current population
		self.sumfit = 0
		self.minfit = 0
		self.avgfit = 0
		self.maxfit = 0
		
		# Generate the initial population
		self.pop = [ Genotype(lchrom, generate) for _ in range(self.popsize) ]
		for gene in self.pop: gene.eval_fitness(obj_func)
		
		# Select arbitrary best gene and update stats to start
		self.best = self.pop[0]
		self.update_statistics()
		return

	# Select a random gene based on fitness
	# Population must be sorted to work effectively
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

	# Select the best gene from a random sample of the population
	def tourn_select(self):
		if self.minimize:
			return min(random.sample(self.pop, TOURN_MAX), key=lambda g: g.fitness)
		else:
			return max(random.sample(self.pop, TOURN_MAX), key=lambda g: g.fitness)
	
	# Evolve the current generation
	def next_generation(self):
		newpop = []
		for _ in range(0, self.popsize, 2): # Do this operation self.popsize/2 times
			# Select two genes to cross
			gene1 = self.tourn_select()
			gene2 = self.tourn_select()
			
			# Cossover our two selected genes
			child1, child2 = self.crossover(gene1, gene2, self.pc)

			newpop.append(child1)
			newpop.append(child2)
		
		# Maintain consistency for odd numbered populations
		if len(newpop) < self.popsize: newpop += [ self.tourn_select() ]
		# Replace old population with new population
		self.pop = newpop

		# Mutate all genes with probability pm
		for i in range(self.popsize):
			# Mutate the gene
			self.pop[i].eval_fitness(self.obj_func)

			# Calculate pmutate based on original pm value; higher chance to mutate worse genes
#			pmutate = self.pm * ( (self.pop[i].fitness - self.minfit) / (self.maxfit - self.minfit) )
#			if self.minimize:
#				pmutate = 1.0 - pmutate
			pmutate = self.pm
			if self.mutate(self.pop[i], pmutate, self.flip_bit):
				self.pop[i].eval_fitness(self.obj_func)
		return

	# Update information on the current generation
	def update_statistics(self):
		# Find total fitness, best gene for this generation
		self.sumfit = sum(g.fitness for g in self.pop)
		self.avgfit = self.sumfit / self.popsize
		
		self.maxfit = max(self.pop, key=lambda x: x.fitness).fitness
		self.minfit = min(self.pop, key=lambda x: x.fitness).fitness
		
		if self.minimize: # Minimize best gene
			self.best = min(self.pop + [self.best], key=lambda x: x.fitness) # Get best gene of generation
		else: # Maximize best gene
			self.best = max(self.pop + [self.best], key=lambda x: x.fitness) # Get best gene of generation
		return

