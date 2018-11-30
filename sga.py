import random
import datetime

MINRNG, MAXRNG = -5.12, 5.12

#Possible options in gen
GENESET = [1, 2, 3, 4, 5]
generation_count = 1

global startTime
startTime = datetime.datetime.now()

# Container class for chromosomes
class Genotype:
    def __init__(self, lchrom, *fitness):
        # Create list of random float values
        temp = []
        while len(temp) < lchrom:
            temp.extend(random.sample(GENESET, len(GENESET)))
        self.chrom = temp  # Encode the values into
        self.fitness = 0
        if fitness:
            self.fitness = int(fitness[0])

    # For printing out gene values
    def __str__(self):
        ret = '{:.3f} | '.format(self.fitness)
        for x in self.chrom:
            ret += '{:.3f} '.format(x)
        return ret

    def set_chrom(self, newchrom):
        self.chrom = newchrom

    def set_fitness(self, newfit):
        self.fitness = newfit


# 0 <= probability <= 1
def flip(probability):
    return random.random() < probability


# Encode a list of floating point values into a chromosome (list of 32 bit int)
def sga_encode(x):
    pass


# Decode a chromosome (list of 32 bit int, but treated like like of booleans) into a list of floats
def sga_decode(chrom):
    pass

# Fitness determination function, checks against target for match
def det_fitness(guess, target):
    total = 0
    for i in range(len(guess)):
        if guess[i] == target[i]:
            total += 1
    return total

# Mutation function, very basic and only changes 1 value in chroma
def mutate(parent, target):
    childgene = parent.chrom[:]

    option1, option2 = random.sample(GENESET, 2)
    index = random.randrange(0, len(parent.chrom))
    if option1 == childgene[index]:
        childgene[index] = option2
    else:
        childgene[index] = option1

    fitness = det_fitness(childgene,target)
    new_gen = Genotype(len(parent.chrom), fitness)
    new_gen.set_chrom(childgene)
    return new_gen


# Determination function to find best strand possible
def get_best(target, max_gen):
    bestparent = Genotype(len(target))
    temp = det_fitness(bestparent.chrom,target)
    bestparent.set_fitness(temp)
    display(bestparent)
    if bestparent.fitness >= len(target):
        return bestparent
    global generation_count
    while generation_count < max_gen:
        random.seed()
        child = mutate(bestparent, target)
        if bestparent.fitness >= child.fitness:
            continue
        generation_count += 1
        display(child)
        if child.fitness >= len(target):
            return child
        bestparent = child

# Output display function for debugging purposes
def display(candidate):
    timeDiff = datetime.datetime.now() - startTime
    print("{} Current Fittness: {} Current Gen: {} Time: {}".format(candidate.chrom, candidate.fitness, generation_count, timeDiff))

# Crossover of two genes
def crossover(chrom1, chrom2, pc):
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
# x => list of decoded values
def himmelblau(x):
    return (x[0] * x[0] + x[1] - 11) ** 2 + (x[0] + x[1] * x[1] - 7) ** 2


if __name__ == '__main__':
    # Input parameters to fiddle with
    maxgen = int(input("Enter maximum # of generations: "))
    lchrom = int(input("Enter chromosome length ----- > "))
    # popsize = int(input("Enter population size ------- > "))
    # pcross = float(input("Enter crossover probability - > "))
    # pmutation = float(input("Enter mutation probability -- > "))

    #Random Target for SGA to determine to
    target = []
    while len(target) < lchrom:
        target.extend(random.sample(GENESET, len(GENESET)))

    print("Target to find: {}".format(target))

    startTime = datetime.datetime.now()
    best = get_best(target, maxgen)
