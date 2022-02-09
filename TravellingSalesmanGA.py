from tkinter import *
from GeneticAlgorithms import *
import math

# Created by @EhsanZaeimzadeh <ehsan@pigmentory.app>
__author__ = "Ehsan Zaeimzadeh"
__email__ = "ehsan@pigmentory.app"
class City:
    def __init__(self, x, y):
        '''
        Representation of a City
        :param x: coordination x
        :param y: coordination y
        :return:
        '''
        self.name = "City [{:0>3} {:0>3}]".format(x, y)
        self.x = x
        self.y = y

    def distanceTo(self, nextCity):
        """
        :param nextCity: next city object in turn
        :return: distance to next city
        """
        xDif = abs(self.x - nextCity.x)
        yDif = abs(self.y - nextCity.y)
        return ( math.sqrt(math.pow(xDif, 2) + math.pow(yDif, 2)) )


class TravellingSalesManGA(GeneticAlgorithms):
    def __init__(self, nrOfCities=10, populationSize=20, mutationProbability=0.05, animationDelay=0.03, realSolution=False):
        """
        :param nrOfCities: Total nr of cities in this game
        :param populationSize: Size of the population
        :param mutationProbability: Probability of mutation to happen
        :param animationDelay: Delay in seconds for animation
        :param realSolution: Whether to calculate the real solution
        """
        GeneticAlgorithms.__init__(self)
        self.nrOfCities = nrOfCities
        self.populationSize = populationSize
        self.mutationProbability = mutationProbability
        self.animationDelay = animationDelay
        self.realSolution = realSolution
        self.mapSize = 620
        self.itterationIndex = 0
        self.setMutationDegree(2)

        self.cities = []
        self.alleles = []
        self.population = []

        self.guiRoot = Tk()
        self.guiRoot.wm_title("Travelling Salesman GA")
        self.guiRoot.resizable(width=False, height=False)
        self.guiDemoF = Frame(self.guiRoot, height=self.mapSize)
        self.guiDemoF.pack()
        self.guiAnimF = Frame(self.guiDemoF, height=self.mapSize, bg="red", relief=GROOVE)
        self.guiLoggF = Frame(self.guiDemoF, height=self.mapSize, bg="red")
        self.guiAnimF.grid(row=0, column=0)
        self.guiLoggF.grid(row=0, column=1)
        self.guiAnimC = Canvas(self.guiAnimF, width=self.mapSize + 20, height=self.mapSize, bg="white", relief=RAISED)
        self.guiAnimC.pack(fill=BOTH)
        self.guiLoggT = Text(self.guiLoggF, width=50, height=40, bg="white")
        self.guiLoggT.pack(fill=BOTH)
        self.guiContF = Frame(self.guiRoot, height=30)
        self.guiContF.pack()
        Button(self.guiContF, text="Evolve1Time!", width=25, relief=GROOVE, bd=4, command=self.evolve1Time).grid(row=0, column=0, padx=2, pady=2)
        Button(self.guiContF, text="Evolve10Time!", width=25, relief=GROOVE, bd=4, command=self.evolve10Time).grid(row=0, column=1, padx=2, pady=2)
        Button(self.guiContF, text="Evolve100Time!", width=25, relief=GROOVE, bd=4, command=self.evolve100Time).grid(row=0, column=2, padx=2, pady=2)
        Button(self.guiContF, text="Evolve1000Time!", width=25, relief=GROOVE, bd=4, command=self.evolve1000Time).grid(row=0, column=3, padx=2, pady=2)
        Button(self.guiContF, text="New Game!", width=25, relief=GROOVE, bd=4, command=self.newGame).grid(row=0, column=5, padx=20, pady=2)
        self.guiCities = []
        self.guiRoutes = []
        self.guiLogges = []
        self.guiText = []

        mainloop()

    def evolve1Time(self):
        self.evolveByNumber(1)
    def evolve10Time(self):
        self.evolveByNumber(10)
    def evolve100Time(self):
        self.evolveByNumber(100)
    def evolve1000Time(self):
        self.evolveByNumber(1000)

    def evolveByNumber(self, nrOfEvolutions):
        self.itterationIndex += nrOfEvolutions
        for i in range(nrOfEvolutions):
            self.population = self.evolve(self.population, Selection_Random, Crossover_OnePoint , Mutation_Swap, self.calculateFitness, True)

        self.guiDrawRoute(self.population.getFittestChromosome().getAlleles(), self.population.getFittestChromosome().getFitness())

    def calculateFitness(self, alleles):
        distance = 0
        for i in range(len(alleles)):
            nextCityIndex = i+1
            if i == len(alleles)-1:
                nextCityIndex = 0
            distance += alleles[i].object.distanceTo(alleles[nextCityIndex].object)
        return float(1/float(distance))

    def newGame(self):
        #Clean up
        self.guiLoggT.delete("1.0", END)
        for guiCity  in self.guiCities:
            self.guiAnimC.delete(guiCity)
        time.sleep(self.animationDelay)
        self.guiCities = []
        self.population = []
        self.itterationIndex = 0

        self.cities = self._createCities()
        self.alleles = self._createAlleles()
        self.population = self._createPopulation()

    def _createPopulation(self):
        chromosomes = []
        fittestChromosome = None
        fittestChromosomeFitness = 0
        for i in range(self.populationSize):
            alleleSample = tuple(random.sample(list(self.alleles), self.nrOfCities))
            fitness = self.calculateFitness(alleleSample)
            chromosome = Chromosome(alleleSample, fitness)
            chromosomes.append(chromosome)
            if fitness > fittestChromosomeFitness:
                fittestChromosomeFitness = fitness
                fittestChromosome = chromosome

        return Population(tuple(chromosomes), fittestChromosome)

    def _createAlleles(self):
        alleles = []
        for city in self.cities:
            alleles.append(Allele(city.name, city))
        return alleles

    def _createCities(self):
        cities = []
        x = list(range(20, self.mapSize-20, 10)); random.shuffle(x)
        y = list(range(20, self.mapSize-20, 10)); random.shuffle(y)
        self.guiLoggT.insert(INSERT, "<<< List of Cities >>>\n")
        for i in range(self.nrOfCities):
            cities.append(City(x[i], y[i]))
            self.guiCities.append(self.guiAnimC.create_rectangle(cities[i].x-2, cities[i].y-2, cities[i].x+2, cities[i].y+2, fill = "black"))
            self.guiLoggT.insert(END, cities[i].name + "\n")
        return cities

    def guiDrawRoute(self, alleles, fitness):
        self.guiRemoveRoute()
        self.guiLoggT.insert(END, "<< Evolution nr: {:<6}   Distance: {}\n".format(self.itterationIndex, str(float(1/fitness))[:8]))
        for i in range(len(alleles)):
            time.sleep(self.animationDelay)
            self.guiAnimC.update()
            self.guiText.append(self.guiAnimC.create_text(alleles[i].object.x-5, alleles[i].object.y-5, text=str(i)))
            nextCityIndex = i+1
            if i == len(alleles)-1:
                nextCityIndex = 0
            self.guiRoutes.append(self.guiAnimC.create_line(alleles[i].object.x,
                                                            alleles[i].object.y,
                                                            alleles[nextCityIndex].object.x,
                                                            alleles[nextCityIndex].object.y,
                                                            fill = "red"))
            time.sleep(self.animationDelay)
            self.guiAnimC.update()

    def guiRemoveRoute(self):
        for r in self.guiRoutes:
            self.guiAnimC.delete(r)
        self.guiRoutes = []
        for t in self.guiText:
            self.guiAnimC.delete(t)
        self.guiText = []



tsGA = TravellingSalesManGA(20)