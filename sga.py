import random

MINRNG, MAXRNG = -5.12, 5.12

# Container class for chromosomes
class Genotype:
	def __init__(self, lchrom):
		# Create list of random float values
		temp = [ random.uniform(MINRNG, MAXRNG) for _ in range(lchrom) ]
		self.chrom = sga_decode( temp ) # Encode the values into 
		self.fitness = 0 # fitness for this gene

	# For printing out gene values
	def __str__(self):
		ret = '{:.3f} | '.format(self.fitness)
		for x in sga_decode(self.chrom):
			ret += '{:.3f} '.format(x)
		return ret

# 0 <= probability <= 1
def flip(probability):
	return random.random() < probability

# Encode a list of floating point values into a chromosome (list of 32 bit int)
def sga_encode(x):
	pass

# Decode a chromosome (list of 32 bit int, but treated like like of booleans) into a list of floats
def sga_decode(chrom):
	pass

# Crossover of two genes
def crossover(chrom1, chrom2, pc):
	pass

# Flip a bit in the chromosome
def mutate(chrom, pm):
	pass
		
# Select a random(ish?) gene
def select(self):
	pass

# Evolve and replace the current population using crossover and mutation
def generation(self):
	pass

# Update statistics for current generation
def statistics(self):
	pass

#  DeJong's Sphere Objective function
# x => list of decoded values
def dejong(x):
	pass

# Rosenbrock's Valley Objective function
# x => list of decoded values
def rosenbrock(x):
	pass

# Himmelblau Objective function
# x => list of decoded values
def himmelblau(x):
	pass

if __name__ == '__main__':
	# Input parameters to fiddle with
	maxgen    =   int(input("Enter maximum # of generations: "))
	popsize   =   int(input("Enter population size ------- > "))
	lchrom    =   int(input("Enter chromosome length ----- > "))
	pcross    = float(input("Enter crossover probability - > "))
	pmutation = float(input("Enter mutation probability -- > "))

	# TODO: Genetic algorithm main
