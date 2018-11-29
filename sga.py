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
    sum = 0
    for i in x:
        sum += i * i
    return sum

# Rosenbrock's Valley Objective function
# x => list of decoded values
def rosenbrock(x):
    sum = 0
    for i in range(len(x) - 1):
        y = x[i]
        y_next = x[i + 1]
        new = 100 * (y_next - y ** 2) ** 2 + (y - 1) ** 2
        sum = sum + new
    return sum

# Himmelblau Objective function NOTE: Apparently this is 2D function (ie.[x,y]) instead of entire list
# x => [x1,x2] length 2 list of decoded values
def himmelblau(x):
    return (x[0] * x[0] + x[1] - 11) ** 2 + (x[0] + x[1] * x[1] - 7) ** 2

if __name__ == '__main__':
	# Input parameters to fiddle with
	maxgen    =   int(input("Enter maximum # of generations: "))
	popsize   =   int(input("Enter population size ------- > "))
	lchrom    =   int(input("Enter chromosome length ----- > "))
	pcross    = float(input("Enter crossover probability - > "))
	pmutation = float(input("Enter mutation probability -- > "))

	# TODO: Genetic algorithm main
