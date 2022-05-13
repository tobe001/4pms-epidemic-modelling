import numpy as np
import matplotlib.pyplot as plt
import random

#Create class for individuals to store their list of neighbours, current state, group, and susceptibility and infectiousness
class Individual:

    def __init__(self, neighbours, state, group, susceptibility, infectiousness):
        self.neighbours = neighbours
        self.state = state
        self.group = group
        self.susceptibility = susceptibility
        self.infectiousness = infectiousness

#Define a function to form a network from a given list of individuals using one of several methods by adding to each individual's list of neighbours (which should initially be empty)
def formNetwork(population, method):

    #Form a ring lattice
    if (method == "ring-lattice"):
        k = 5 #Parameter of method
        #For each individual i (consider them placed in a ring):
        for i in range(0, len(population)):
            #Connect i to the k following individuals, looping around to individual 0 when i is near the end of the population
            for j in range(1, k + 1):
                (population[i].neighbours).append(population[(i + j) % len(population)])
                (population[(i + j) % len(population)].neighbours).append(population[i])

    #Form a binomial random network
    elif (method == "binomial"):
        p = 0.01 #Parameter of method
        #For each pair of individuals i and j such that j > i (this ensures that each pair is considered only once):
        for i in range(0, len(population)):
            for j in range(i + 1, len(population)):
                #Link i and j with probability p
                rand = np.random.uniform(0, 1)
                if (rand <= p):
                    (population[i].neighbours).append(population[j])
                    (population[j].neighbours).append(population[i])

    #Form a small-world network
    elif (method == "small-world"):
        k = 5 #Parameter of method
        p = 0.1 #Parameter of method
        #For each individual i (consider them placed in a ring):
        for i in range(0, len(population)):
            #Form a booleam list of length k which will determine which of the edges to the k individuals following i will be rewired
            rewire = np.zeros(k)
            for j in range(0, k):
                #Mark each edge for rewiring with probability p
                rand = np.random.uniform(0, 1)
                if (rand <= p):
                    rewire[j] = True
                else:
                    rewire[j] = False
            #For j such that we do not rewire the edge from i to the jth individual ahead of it:
            for j in range(1, k + 1):
                if not(rewire[j - 1]):
                    #Link i to the jth individual ahead of it, looping around to indivdual 0 when i is near the end of the population
                    (population[i].neighbours).append(population[(i + j) % len(population)])
                    (population[(i + j) % len(population)].neighbours).append(population[i])
            #For j such that we do rewire the edge from i to the jth individual ahead of it:
            for j in range(1, k+1):
                if rewire[j - 1]:
                    #Choose an individual uniformlly at random, such that the chosen individual is not i or the jth individual ahead of i, and is not already linked to i
                    chosen = np.random.randint(0, len(population))
                    while ((chosen == i) or (chosen == ((i + j) % len(population))) or (population[chosen] in population[i].neighbours)):
                        chosen = np.random.randint(0, len(population))
                    #Link i and the chosen individual
                    (population[i].neighbours).append(population[chosen])
                    (population[chosen].neighbours).append(population[i])

    #Form a scale-free network
    elif (method == "scale-free"):
        m = 5 #Parameter of method
        #Form fully connected network from the first m individuals
        for i in range(0, m):
            for j in range(i+1, m):
                (population[i].neighbours).append(population[j])
                (population[j].neighbours).append(population[i])
        #To add each remaining individual:
        for i in range(m, len(population)):
            degrees = np.zeros(i)
            probs = np.zeros(i)
            cumProbs = np.zeros(i)
            #Calculate degree of each individual currently in the network and the sum of all degrees
            for j in range(0, i):
                degrees[j] = len(population[j].neighbours)
            degSum = sum(degrees)
            #Calculate the probability of individual i connecting to each individual based on the formula in the BA model, and the cumulative probabilities
            for j in range(0, i):
                probs[j] = degrees[j]/degSum
                cumProbs[j] = sum(probs[0 : j + 1])
            #Select the m existing individuals to which individual i will initially connect
            newNeighbours = []
            for k in range(0, m):
                selected = False
                while not selected:
                    rand = np.random.uniform(0, 1)
                    #For each existing individual j, if j is the first individual for which our randomly generated number is less than the associated cumulative probability, then select individual j (provided j has not already been selected)
                    #This is because the gap between cumProbs[j-1] and cumProbs[j] will be the probability of connecting to individual j
                    for j in range(0, i):
                        if (rand < cumProbs[j]):
                            if not(j in newNeighbours):
                                newNeighbours.append(j)
                                selected = True
                            break
            #Connect individual i to the m selected individuals
            for k in range(0, m):
                (population[i].neighbours).append(population[newNeighbours[k]])
                (population[newNeighbours[k]].neighbours).append(population[i])
    
    return population  

#Define a function to return a list of inidces of individuals in the population from a given subset sorted by one of several measures
def sortedSubsetofPopulation(population, subset, measure):

    if (len(subset) > 0):
        sortedSubset = [subset[0]]
    else:
        return []
    
    #Sort individuals randomly
    if (measure == "random"):
        #For each individual in the subset:
        for index in range(1, len(subset)):
            label = subset[index]
            #Generate random index in list of sorted individuals (including end of list) and insert current individual at that index
            i = np.random.randint(0, len(sortedSubset) + 1)
            sortedSubset.insert(i, label)

    #Sort individuals by degree
    elif (measure == "degree"):
        #For each individual in the subset:
        for index in range(1, len(subset)):
            label = subset[index]
            currentIndividual = population[label]
            currentDegree = len(currentIndividual.neighbours)
            #Compate current individual's degree with degree of each individual that has already been sorted
            for j in range(0 , len(sortedSubset)):
                comparisonIndividual = population[sortedSubset[j]]
                comparisonDegree = len(comparisonIndividual.neighbours)
                #If current individual's degree is greater than or equal to that of the individual it is being compared against, insert the former into the list ahead of the latter
                if (currentDegree >= comparisonDegree):
                    sortedSubset.insert(j, label)
                    break
                #If the end of the list has been reached and the current individual has not been inserted, add them to the end of the list
                if (j == (len(sortedSubset) - 1)):
                    sortedSubset.append(label)

    #Sort individuals by susceptibility
    elif (measure == "susceptibility"):
        #For each individual in the subset:
        for index in range(1, len(subset)):
            label = subset[index]
            currentIndividual = population[label]
            #Compare current individual's susceptibility with susceptibility of each individual that has already been sorted
            for j in range(0 , len(sortedSubset)):
                comparisonIndividual = population[sortedSubset[j]]
                #If current individual's susceptibility is greater than or equal to that of the individual it is being compared against, insert the former into the list ahead of the latter
                if (currentIndividual.susceptibility >= comparisonIndividual.susceptibility):
                    sortedSubset.insert(j, label)
                    break
                #If the end of the list has been reached and the current individual has not been inserted, add them to the end of the list
                if (j == (len(sortedSubset) - 1)):
                    sortedSubset.append(label)

    #Sort individuals by infectiousness
    elif (measure == "infectiousness"):
        #For each individual in the subset:
        for index in range(1, len(subset)):
            label = subset[index]
            currentIndividual = population[label]
            #Compare current individual's infectiousness with infectiousness of each individual that has already been sorted
            for j in range(0 , len(sortedSubset)):
                comparisonIndividual = population[sortedSubset[j]]
                #If current individual's infectiousness is greater than or equal to that of the individual it is being compared against, insert the former into the list ahead of the latter
                if (currentIndividual.infectiousness >= comparisonIndividual.infectiousness):
                    sortedSubset.insert(j, label)
                    break
                #If the end of the list has been reached and the current individual has not been inserted, add them to the end of the list
                if (j == (len(sortedSubset) - 1)):
                    sortedSubset.append(label)

    #Sort individuals by group only (with group B taking priority over A)
    elif (measure == "group"):
        #Create seperate lists for individuals in groups A and B
        groupA = []
        groupB = []
        #For each individual in the subset:
        for index in range(0, len(subset)):
            label = subset[index]
            currentIndividual = population[label]
            #Add individual to list of corresponding group
            if (currentIndividual.group == "A"):
                groupA.append(label)
            elif (currentIndividual.group == "B"):
                groupB.append(label)
        #Concatenate lists, with the list for group B coming first
        sortedSubset = groupB + groupA
    
    return sortedSubset

#Set parameters
N = 1000
I0 = 1
R0 = 0
p = 0.7
gamma = 0.3
mu = 0.1
subsetSize = N
maxT = 120
numSims = 100

vs = np.linspace(0, 1000, 101)

#Set random seeds for reproducability
np.random.seed(1)
random.seed(1)

#Initialise individuals
population = []
for i in range(0, N):
    #Assign each individual to group A or B with equal probability and generate their susceptibility and infectiousness accordingly
    rand = np.random.uniform(0, 1)
    if (rand <= 0.5):
        group = "A"
        sus = np.random.beta(6, 6)
        inf = np.random.beta(4, 2)
    else:
        group = "B"
        sus = np.random.beta(10, 2)
        inf = np.random.beta(4, 2)
    population.append(Individual([], "Uninitialised", group, sus, inf))
#Form network
population = formNetwork(population, "small-world")

#Initialise 3D list to store fatality numbers for all runs for all values of v
listofListsofDs = []

#For each value of v:
for v in vs:

    #Initialise 2D list to hold fatality numbers for all runs
    listofDs = []

    #Run simulation multiple times
    for n in range(0, numSims):
    
        #Initialise fatality numbers
        D = np.zeros(maxT + 1)
        D[0] = 0

        #Initialise individual states
        #Make first S0 individuals susceptible
        for i in range(0, N - I0 - R0):
            population[i].state = "S"
        #Make next I0 individuals infectious
        for i in range(N - I0 - R0, N - R0):
            population[i].state = "I"
        #Make next R0 individuals recovered
        for i in range(N - R0, N):
            population[i].state = "R"

        #Generate random subset of individuals that can be vaccinated
        potentialVaccineTakers = random.sample(range(0, N), subsetSize)
        #Sort individuals in subset by chosen measure of importance
        sortedPotentialVaccineTakers = sortedSubsetofPopulation(population, potentialVaccineTakers, "random")

        #Vaccinate v most important susceptible individuals from the subset of potential vaccine takers
        vaccinesAvailable = v
        index = 0
        #While we have not yet used all the vaccines and we have not yet considered all individuals in the subset of potential takers:
        while ((vaccinesAvailable > 0) and (index < subsetSize)):
            #Consider next individual in the sorted subset
            currentIndividual = population[sortedPotentialVaccineTakers[index]]
            #If current individual is susceptible, vaccinate them (make them recovered and reduce the number of available vaccines by one)
            if (currentIndividual.state == "S"):
                currentIndividual.state = "R"
                vaccinesAvailable -= 1
            index += 1

        #Run simulation
        for t in range(1, maxT + 1):
            #Create list to store states of individuals for next time step
            nextStates = [""] * N
            #For each individual:
            for i in range(0, N):
                currentIndividual = population[i]
            
                #Determine next state for susceptible individuals
                if (currentIndividual.state == "S"):
                    infected = False
                    #For each neighbour of the current individual:
                    for j in range(0, len(currentIndividual.neighbours)):
                        currentNeighbour = (currentIndividual.neighbours)[j]
                        #If the neighbour is infectious, they infect the current individual with probability p * current individual's suceptibility * neighbour's infectiousness
                        if (currentNeighbour.state == "I"):
                            rand = np.random.uniform(0, 1)
                            if (rand <= p * currentIndividual.susceptibility * currentNeighbour.infectiousness):
                                infected = True
                                #If current individual is infected by one neighbour, there is no need to check if they are infected by any others
                                break
                    if infected:
                        nextStates[i] = "I"
                    else:
                        nextStates[i] = "S"
                
                #Determine next state for infectious individuals
                elif (currentIndividual.state == "I"):
                    #Recovers with probability gamma
                    rand = np.random.uniform(0, 1)
                    if (rand <= gamma):
                        nextStates[i] = "R"
                    #Dies with probability mu
                    elif (rand <= gamma + mu):
                        nextStates[i] = "D"
                    #Otherwise, remains infectious
                    else:
                        nextStates[i] = "I"

                #Determine next state for recovered individuals
                elif (currentIndividual.state == "R"):
                    #Always remains recovered
                    nextStates[i] = "R"

                #Determine next state for deceased individuals
                elif (currentIndividual.state == "D"):
                    #Always remains deceased
                    nextStates[i] = "D"

            #Update states and count number of dead
            deceasedCount = 0
            for i in range(0, N):
                population[i].state = nextStates[i]
                if (population[i].state == "D"):
                    deceasedCount += 1
            D[t] = deceasedCount

        #Check to see if infection died out in this run - if not, we must simulate for greater period of time
        for i in range(0, N):
            if (population[i].state == "I"):
                print("ERROR: Infection did not die out when v = " + str(v) + ". Try increasing maxT.")

        #Store fatality numbers for this run
        listofDs.append(D)
        print("Completed simulation " + str(n) + " for v = " + str(v) + ".")

    #Store list of fatality numbers for this v
    listofListsofDs.append(listofDs)

#Calculate average final fatality numbers and standard deviation for each v
averageFinalDs = np.zeros(len(vs))
stdDevFinalDs = np.zeros(len(vs))
for k in range(0, len(vs)):
    averageFinalDs[k] = np.mean([listofListsofDs[k][n][maxT] for n in range(0, numSims)])
    stdDevFinalDs[k] = np.std([listofListsofDs[k][n][maxT] for n in range(0, numSims)])

plt.rcParams.update({'font.size': 14})
plt.figure()
#Plot average final fatality numbers against v
plt.plot(vs, averageFinalDs, 'k-', linewidth = 2.5)
plt.fill_between(vs, averageFinalDs - stdDevFinalDs, averageFinalDs + stdDevFinalDs, color = 'blue', alpha = 0.3)
plt.xlim([0, 1000])
plt.ylim([0, 300])
plt.xlabel('v')
plt.ylabel('Final number of deceased individuals')

plt.show()
