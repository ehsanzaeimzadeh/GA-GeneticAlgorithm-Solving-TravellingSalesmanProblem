import random, datetime, time
# Created by @EhsanZaeimzadeh <ehsan@pigmentory.app>
__author__ = "Ehsan Zaeimzadeh"
__email__ = "ehsan@pigmentory.app"


Selection_RouletteWheel = "RouletteWheel"
Selection_StochasticUniversal = "StochasticUniversal"
Selection_Tournament = "Tournament"
Selection_Random = "Random"
Crossover_OnePoint = "OnePoint"
Crossover_MultiPoint = "MultiPoint"
Crossover_Uniform = "UniForm"
Crossover_OddEven = "OddEven"
Mutation_Swap = "Swap"
Mutation_Random = "Random"
Mutation_Scramble = "Scramble"
Mutation_Invert = "Invert"

#-----------------------------------------------------------------------------------------------------------------------
# Allele
#-----------------------------------------------------------------------------------------------------------------------
class Allele(object):
    def __init__(self, id,  object, fitness = 0.0):
        self.id = str(id)
        self.object = object
        self.fitness = fitness

    def __str__(self):
        return self.object.__str__()
#-----------------------------------------------------------------------------------------------------------------------
# Chromosome
#-----------------------------------------------------------------------------------------------------------------------
class Chromosome(object):
    def __init__(self, alleles, fitness = 0):
        self._alleles = tuple(alleles)
        self._fitness = fitness
        self._alleleIdMAPindex = {}
        i = 0
        for allele in self._alleles:
            self._alleleIdMAPindex[allele.id] = i
            i+=1

    def getAlleles(self): return list(self._alleles)
    def getAllele(self, gene): return self._alleles[gene]
    def getAlleleById(self, id):
        try:
            a = self._alleles[self._alleleIdMAPindex[id]]
        except:
            a = None
        finally:
            return a
    def getSize(self): return len(self._alleles)
    def setFitness(self, fitness): self._fitness = fitness
    def getFitness(self): return self._fitness

    def display(self, seperator = "-", doPrint = True):
        text =  "Alleles: "
        for allele in self._alleles:
            text += allele.id + seperator
        text += " Fitness: " + str(self._fitness)
        if doPrint: print(text)
        return text
#-----------------------------------------------------------------------------------------------------------------------
# Population
#-----------------------------------------------------------------------------------------------------------------------
class Population(object):
    def __init__(self, chromosomes, fittestChromosome = None):
        self._chromosomes = tuple(chromosomes)
        self._fittest = fittestChromosome

    def getChromosomes(self): return list(self._chromosomes)
    def getChromosome(self, index): return self._chromosomes[index]
    def getSize(self): return len(self._chromosomes)
    def setFittestChromosome(self, chromosome): self._fittest = chromosome
    def getFittestChromosome(self): return self._fittest
#-----------------------------------------------------------------------------------------------------------------------
# GeneticAlgorithms
#-----------------------------------------------------------------------------------------------------------------------
class GeneticAlgorithms(object):
    def __init__(self, mutationProbability = 0.05, selectionSize = 5):
        """

        :param mutationProbability:
        :param selectionSize:
        """
        self._mutationProbability = mutationProbability
        self._mutationDegree = 5
        self._noDuplicates = True
        self._log = False
        self._selectionSize = selectionSize
        self._logFile = None
        self._alleles = None
        self.selectionMethods = {"RouletteWheel":       "chromosomeSelection_RouletteWheel",
                                 "StochasticUniversal": "chromosomeSelection_StochasticUniversal",
                                 "Tournament":          "chromosomeSelection_Tournament",
                                 "Random":              "chromosomeSelection_Random"}
        self.crossoverMethods = {"OnePoint":            "chromosomeCrossover_OnePoint",
                                 "MultiPoint":          "chromosomeCrossover_MultiPoint",
                                 "UniForm":             "chromosomeCrossover_Uniform",
                                 "OddEven":             "chromosomeCrossover_OddEven"}
        self.mutationMethods  = {"Swap":                "chromosomeMutation_Swap",
                                 "Random":              "chromosomeMutation_RandomResetting",
                                 "Scramble":            "chromosomeMutation_Scramble",
                                 "Invert":              "chromosomeMutation_Invert"}

    def getMutationProbability(self): return self._mutationProbability
    def setMutationProbability(self, mutationProbability): self._mutationProbability = mutationProbability
    def getSelectionSize(self): return self._selectionSize
    def setSelectionSize(self, selectionSize): self._selectionSize = selectionSize
    def setAlleles(self, alleles): self._alleles = alleles
    def setMutationDegree(self, mutationDegree): self._mutationDegree = mutationDegree
    def setNoDuplicates(self, noDuplicates): self._noDuplicates = bool(noDuplicates)

    def openLog(self):
        self._log = True
        dateTag = str(datetime.datetime.now()).replace(":", "").replace(" ", "").replace("-", "")[:14]
        self._logFile = open("LogFile_{}.txt".format(dateTag), "a")

    def closeLog(self):
        self._logFile.close()

    def _logAlleles(self, alleles, fitness, indent = 0, header = ""):
        if header != "":
            self._logFile.write("{}\n".format(header))
        self._logFile.write(" "*indent + "Alleles: {} Fitness: {}\n".format("".join([allele.id + "-" for allele in alleles]), fitness))

    def _logText(self, text):
        self._logFile.write(text + "\n")

    def getFittestChromosome(self, chromosomes):
        fitness0 = 0
        fittest = None
        for chromosome in chromosomes:
            fitness1 = chromosome.getFitness()
            if fitness1 >= fitness0:
                fitness0 = fitness1
                fittest = chromosome
        return fittest

    def sortChromosomesByFitness(self, chromosomes):
        sortedChromosomes = []
        chromosomesDict = {chromosome.getFitness():chromosome for chromosome in chromosomes}
        for k in reversed(sorted(chromosomesDict.keys())):
            sortedChromosomes.append(chromosomesDict[k])
        return sortedChromosomes

    def chromosomeSelection_RouletteWheel(self, population):
        selectedChromosome = Chromosome([Allele("None", None), Allele("None", None)], 0)
        fitnessSum0 = 0
        for chromosome in population.getChromosomes():
            fitnessSum0 += chromosome.getFitness()

        for i in range(self._selectionSize):
            threshold = random.uniform(0, fitnessSum0)
            fitnessSum1 = 0
            for chromosome in population.getChromosomes():
                fitnessSum1 += chromosome.getFitness()
                if fitnessSum1 >= threshold:
                    if chromosome.getFitness() >= selectedChromosome.getFitness():
                        selectedChromosome = chromosome
        return selectedChromosome

    def chromosomeSelection_StochasticUniversal(self, population):
        a = random.randrange(0, population.getSize() - self._selectionSize)
        return self.getFittestChromosome(population.getChromosomes()[a:a+self._selectionSize])

    def chromosomeSelection_Tournament(self, population):
        selection = random.sample(population.getChromosomes(), self._selectionSize)
        return self.getFittestChromosome(selection)

    def chromosomeSelection_Random(self, population):
        selected = self.getFittestChromosome(random.sample(population.getChromosomes(), 1))
        return selected

    def chromosomeCrossover_OnePoint(self, chromosome1, chromosome2):
        """
        :param chromosome1:
        :param chromosome2:
        :return:
        """
        offspringAlleles = []
        threshold = random.randint(0, chromosome1.getSize()-1)
        if self._noDuplicates:
            offspringAlleles = chromosome1.getAlleles()[:threshold]
            for allele in chromosome2.getAlleles():
                if allele not in offspringAlleles:
                    offspringAlleles.append(allele)
                    if len(offspringAlleles) == chromosome1.getSize():
                        break
        else:
            offspringAlleles = chromosome1.getAlleles()[:threshold] + chromosome2.getAlleles()[threshold:]

        return offspringAlleles

    def chromosomeCrossover_MultiPoint(self, chromosome1, chromosome2):
        offspringAlleles = []
        a, b = sorted(random.sample(range(chromosome1.getSize()), 2))
        if self._noDuplicates:
            offspringAlleles = chromosome1.getAlleles()[a:b]
            for allele in chromosome2.getAlleles():
                if allele not in offspringAlleles:
                    offspringAlleles.append(allele)
                    if len(offspringAlleles) == chromosome1.getSize():
                        break
            if len(offspringAlleles) != chromosome1.getSize():
                for allele in chromosome1.getAlleles():
                    if allele not in offspringAlleles:
                        offspringAlleles.append(allele)
                        if len(offspringAlleles) == chromosome1.getSize():
                            break
        else:
            offspringAlleles = chromosome1.getAlleles()[:a] + chromosome2.getAlleles()[a:b] + chromosome1.getAlleles()[b:]

        return offspringAlleles

    def chromosomeCrossover_Uniform(self, chromosome1, chromosome2):
        alleleIds = [a.id for a in chromosome1.getAlleles()] + [a.id for a in chromosome2.getAlleles()]
        idSamples = random.sample(alleleIds, chromosome1.getSize())
        offspringAlleles = []
        for i in range(chromosome1.getSize()):
            if chromosome1.getAlleleById(idSamples[i]) == None:
                offspringAlleles.append(chromosome2.getAlleleById(idSamples[i]))
            else:
                offspringAlleles.append(chromosome1.getAlleleById(idSamples[i]))
        #offspringAlleles = random.sample(chromosome1.getAlleles() + chromosome2.getAlleles(), chromosome1.getSize())
        print ("\nChromosome 1")
        chromosome1.display()

        print ("\nChromosome 2")
        chromosome2.display()
        print ("\nOffspring")
        for allele in offspringAlleles:
            print (allele.id,)
        print ("\n\n")
        return offspringAlleles

    def chromosomeCrossover_OddEven(self, chromosome1, chromosome2):
        offspringAlleles = []
        if self._noDuplicates:
            for i in range(chromosome1.getSize()):
                if i%2 == 0:
                    if chromosome1.getAllele(i) not in offspringAlleles:
                        offspringAlleles.append(chromosome1.getAllele(i))
                else:
                    if chromosome2.getAllele(i) not in offspringAlleles:
                        offspringAlleles.append(chromosome2.getAllele(i))
            if len(offspringAlleles)!= chromosome1.getSize():
                for i in range(chromosome1.getSize()):
                    if i % 2 == 0:
                        if chromosome2.getAllele(i) not in offspringAlleles:
                            offspringAlleles.append(chromosome2.getAllele(i))
                    else:
                        if chromosome1.getAllele(i) not in offspringAlleles:
                            offspringAlleles.append(chromosome1.getAllele(i))
                    if len(offspringAlleles) == chromosome1.getSize():
                        break
        else:
            offspringAlleles = [chromosome1.geAllele(i) if i%2==0 else chromosome2.getAllele(i) for i in range(chromosome1.getSize())]
        return offspringAlleles

    def chromosomeMutation_Swap(self, chromosomeAlleles):
        alleles = chromosomeAlleles
        if random.random() < self._mutationProbability:
            for i in range(self._mutationDegree):
                a, b = random.sample(range(len(alleles)), 2)
                alleles[a], alleles[b] = alleles[b], alleles[a]
        return alleles

    def chromosomeMutation_RandomResetting(self, chromosomeAlleles):
        alleles = chromosomeAlleles
        if random.random() < self._mutationProbability:
            for i in range(self._mutationDegree):
                alleles[random.randrange(0, len(alleles))] = random.sample(set(self._alleles) - set(alleles), 1)[0]
        return alleles

    def chromosomeMutation_Scramble(self, chromosomeAlleles):
        alleles = chromosomeAlleles
        if random.random() < self._mutationProbability:
            position = random.randrange(0, len(alleles)-self._mutationDegree)
            subset = alleles[position:position+self._mutationDegree]
            random.shuffle(subset)
            alleles[position:position+self._mutationDegree] = subset
        return alleles

    def chromosomeMutation_Invert(self, chromosomeAlleles):
        alleles = chromosomeAlleles
        if random.random() < self._mutationProbability:
            position = random.randrange(0, len(alleles) - self._mutationDegree)
            subset = alleles[position:position + self._mutationDegree]
            subset.reverse()
            alleles[position:position + self._mutationDegree] = subset
        return alleles

    def evolve(self, population, selectionMethod, crossoverMethod, mutationMethod, fitnessFunction, elitism = False):
        t1 = time.time()
        if self._log:
            self._logText("\n" + "-"*100+"\nIteration Start\nSelection: {}    Crossover: {}    Mutation: {}\n".format(selectionMethod, crossoverMethod, mutationMethod)+"-"*100)
        chromosomes = [population.getFittestChromosome()]
        fittestChromosome = None
        fittestChromosomeFitness = 0
        startIndex = 0
        if elitism:
            startIndex = 1
            fittestChromosome = population.getFittestChromosome()
            fittestChromosomeFitness = population.getFittestChromosome().getFitness()
        for i in range(startIndex, population.getSize()):
            motherChromosome = getattr(self, self.selectionMethods[selectionMethod])(population)
            if self._log:
                self._logAlleles(motherChromosome.getAlleles(), motherChromosome.getFitness(), 3, "Mother/Father/Child/Mutated")
            fatherChromosome = getattr(self, self.selectionMethods[selectionMethod])(population)
            if self._log:
                self._logAlleles(fatherChromosome.getAlleles(), fatherChromosome.getFitness(), 3, "")
            childAlleles = getattr(self, self.crossoverMethods[crossoverMethod])(motherChromosome, fatherChromosome)
            if self._log:
                self._logAlleles(childAlleles, "N/A", 3, "")
            mutatedAlleles = getattr(self, self.mutationMethods[mutationMethod])(childAlleles)
            chromosomeFitness = fitnessFunction(mutatedAlleles)
            if self._log:
                self._logAlleles(mutatedAlleles, chromosomeFitness, 3, "")
            chromosome = Chromosome(mutatedAlleles, chromosomeFitness)
            chromosomes.append(chromosome)
            if chromosomeFitness >= fittestChromosomeFitness:
                fittestChromosomeFitness = chromosomeFitness
                fittestChromosome = chromosome
        if self._log:
            self._logAlleles(population.getFittestChromosome().getAlleles(), population.getFittestChromosome().getFitness(), 3, "Fittest Chromosome (before):")
            self._logAlleles(fittestChromosome.getAlleles(), fittestChromosome.getFitness(), 3, "Fittest Chromosome (after):")
            self._logText("-" * 100 + "\nIteration Stop    -    Elapsed Time: {}\n".format(time.time() - t1) + "-" * 100 + "\n")

        return Population(tuple(chromosomes), fittestChromosome)

#-----------------------------------------------------------------------------------------------------------------------