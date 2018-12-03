import random
import matplotlib.pyplot as plt
from copy import deepcopy

from objfunc import *
from sga import SGA, flip

# Flip a random bit in a given chromosome
def flip_bit_bcd(chrom):
	i = random.randint(0, len(chrom)-1)
	x = chrom[i]
	x = int(x * (10**PRECISION)) # Convert to integer (keeping precision)
	sign = -1 if x < 0 else 1 # Get sign of number
	
	# Give sign bit a 1/4 chance compared to full digits (BCD)
	k = random.randint(-1, (NUM_DIGITS-1) * 4)

	# Flip sign bit and return
	if k < 0: 
		chrom[i] = -x / (10**PRECISION)
		return
	else:
		k //= 4
	
	x = abs(x)
	digit = (x // (10**k)) % 10 # Get random digit k
	x -= digit * (10**k)        # Clear digit k

	# Select only values which will result in a valid output
	exclusions = [ j for j in range(0, 9+1) if not (MINRNG <= (x + j*(10**k)) / (10**PRECISION) <= MAXRNG) ]

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
	if len(opts) == 0: 
		return # No options, no change

	# Set digit k to new digit
	ret = x + random.choice( opts ) * (10**k)
	ret = sign * (ret / (10**PRECISION))
	chrom[i] = ret
	return

def flip_bit_plus_minus(chrom):
	i = random.randint(0, len(chrom)-1)
	chrom[i] = - chrom[i]
	return
	
def flip_bit_nqueens(chrom):
	n = int(len(chrom)**0.5)
	i = random.randint(0, n-1) * n
	
	q = [ j for j in range(i, i+n) if chrom[j] ]
	qnot = [ j for j in range(i, i+n) if not chrom[j] ]
	
	i1 = random.choice(q)
	i2 = random.choice(qnot)
	
	chrom[i1] = 0 # Remove queen from i1
	chrom[i2] = 1 # Place queen at i2
	return

def crossover_nqueens(p1, p2, pc):
	c1, c2 = deepcopy(p1), deepcopy(p2)
	if flip(pc):
		n = int(len(c1.chrom)**0.5)
		i = random.randint(0, n-1) * n
		c1.chrom[:i], c2.chrom[i:] = c2.chrom[:i], c1.chrom[i:]
	return c1, c2

if __name__ == '__main__':
	maxgen    =   int(input("Enter maximum # of generations: "))
	popsize   =   int(input("Enter population size ------- > "))
	pcross    = float(input("Enter crossover probability - > "))
	pmutation = float(input("Enter mutation probability -- > "))
	minimize  =       input("Minimize? (Y/N) ------------- > ") # If user answers No, maximize instead
	minimize  = (minimize[0].lower() == 'y') # Get minimize as boolean

	# GLOBAL VARIABLES
	MINRNG, MAXRNG = -5.12, 5.12 # Min and max for floating point numbers
	NUM_DIGITS = 3 # How many significant digits for each value
	PRECISION = NUM_DIGITS-1 # How many decimal places for floating points (must be <= NUM_DIGITS)
	CONVERGENCE = 200 # If no improvement after so many generations, exit

	########################### MAIN CODE : RUN TEST CASES FROM HERE ############################
#	lchrom = 27 # for objfuncXX, set lchrom = XX
#	objfunc = objfunc27
#	generate = lambda x: [random.choice( [-1, 1] ) for _ in range(x)] # lambda for selecting random chromosome
#	sga = SGA(popsize, lchrom, pmutation, pcross, objfunc, generate, flip_bit_plus_minus, minimize)

	lchrom = 30
	objfunc = dejong
	generate = lambda x: [random.uniform( MINRNG, MAXRNG ) for _ in range(x)]
	sga = SGA(popsize, lchrom, pmutation, pcross, objfunc, generate, flip_bit_bcd, minimize)

	'''
	def gen_queens(n): # Function to generate n queens
		queens = []
		for i in range(n):
			t = [0] * n
			t[random.randint(0, n-1)] = 1
			queens += t
		return queens

	lchrom = int(input("Enter # of queens: "))
	objfunc = nqueens
	generate = gen_queens
	sga = SGA(popsize, lchrom, pmutation, pcross, objfunc, generate, flip_bit_nqueens, minimize, cross=crossover_nqueens)
	'''
	print() # Add extra line for prettyness
	#############################################################################################

	# The main body of the SGA, evolve to next generation and update statistics until max no. of generations have been parsed
	y = []
	lastfit, converge = sga.best.fitness, 0
	for gen in range(maxgen):
		sga.next_generation()
		sga.update_statistics()
		
		y.append(sga.best.fitness)
		if sga.best.fitness - lastfit == 0:
			converge += 1
			if converge == CONVERGENCE:
				print("Convergence reached at generation {:d}".format(gen))
				break
		else:
			converge = 0
		lastfit = sga.best.fitness
		
	print("Best Gene:", str(sga.best))
	
	x = list(range(gen+1))
	plt.plot(x, y, label='best')
	plt.xlabel('Generation #')
	plt.ylabel('Best, Avg Fitness')
	plt.show()

