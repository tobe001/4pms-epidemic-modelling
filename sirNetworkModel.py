import numpy as np
import matplotlib.pyplot as plt

#Create class for individuals to store their list of neighbours and current state
class Individual:

    def __init__(self, neighbours, state):
        self.neighbours = neighbours
        self.state = state

#Define function to form network from a given list of individuals by one of several methods
def formNetwork(population, method):

    #Form a rectangular lattice
    if (method == "rectangular-lattice"):
        ncol = 50 #Change as parameter
        #For each pair of individuals i and j such that j > i (consider them placed in a rectangular grid from left to right, then top to bottom):
        for i in range(0, len(population)):
            for j in range(i + 1, len(population)):
                #If j is 1 ahead of i and i is not on the end of a row, link i and j
                if (((i % ncol) != (ncol - 1)) and ((j - i) == 1)):
                    (population[i].neighbours).append(population[j])
                    (population[j].neighbours).append(population[i])
                #If j is one row ahead of i, link i and j
                if ((j - i) == ncol):
                    (population[i].neighbours).append(population[j])
                    (population[j].neighbours).append(population[i])

    #Form a ring lattice
    elif (method == "ring-lattice"):
        k = 5 #Change as parameter
        #For each individual i (consider them placed in a ring):
        for i in range(0, len(population)):
            #Connect i to the k following individuals, looping around to individual 0 when i is near the end of the population
            for j in range(1, k+1):
                (population[i].neighbours).append(population[(i+j) % len(population)])
                (population[(i+j) % len(population)].neighbours).append(population[i])

    #Form a network with random connections
    elif (method == "random"):
        p = 0.01 #Change as parameter
        #For each pair of individuals i and j such that j > i:
        for i in range(0, len(population)):
            for j in range(i + 1, len(population)):
                #Link i and j with probability p
                rand = np.random.uniform(0, 1)
                if (rand <= p):
                    (population[i].neighbours).append(population[j])
                    (population[j].neighbours).append(population[i])

    #Form a small-world network
    elif (method == "small-world"):
        k = 5 #Change as parameter
        p = 0.1 # Change as parameter
        #For each individual i (consider them placed in a ring):
        for i in range(0, len(population)):
            #Form a list of length k determining which of the edges to the k individuals following i will be rewired
            rewire = np.zeros(k)
            for j in range(0, k):
                #Rewire each edge with probability p
                rand = np.random.uniform(0, 1)
                if (rand <= p):
                    rewire[j] = True
                else:
                    rewire[j] = False
            #For j such that we do not rewire the edge from i to the jth individual ahead of it:
            for j in range(1, k+1):
                if not(rewire[j-1]):
                    #Link i to the jth individual ahead of it, looping around to indivdual 0 when i is near the end of the population
                    (population[i].neighbours).append(population[(i+j) % len(population)])
                    (population[(i+j) % len(population)].neighbours).append(population[i])
            #For j such that we do rewire the edge from i to the jth individual ahead of it:
            for j in range(1, k+1):
                if rewire[j-1]:
                    #Choose an individual uniformlly at random, such that the chosen individual is not i or the jth individual ahead of i, and is not already linked to i
                    chosen = np.random.randint(0, len(population))
                    while ((chosen == i) or (chosen == ((i+j) % len(population))) or (population[chosen] in population[i].neighbours)):
                        chosen = np.random.randint(0, len(population))
                    #Link i and the chosen individual
                    (population[i].neighbours).append(population[chosen])
                    (population[chosen].neighbours).append(population[i])

    #Form a scale-free network
    elif (method == "scale-free"):
        m = 5 #Change as parameter
        #Form fully connected network of the first m individuals
        for i in range(0, m):
            for j in range(i+1, m):
                (population[i].neighbours).append(population[j])
                (population[j].neighbours).append(population[i])
        #To add each remaining individual:
        for i in range(m, len(population)):
            degrees = np.zeros(i)
            probs = np.zeros(i)
            cumProbs = np.zeros(i)
            #Calculate degrees of each individual currently in the network and the sum of all degrees
            for j in range(0, i):
                degrees[j] = len(population[j].neighbours)
            degSum = sum(degrees)
            #Calculate the probability of individual i connecting to each individual based on the formula in the BA model, and the cumulative probabilities
            for j in range(0, i):
                probs[j] = degrees[j]/degSum
                cumProbs[j] = sum(probs[0:j+1])
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

#Define parameters
N = 1000
I0 = 1
R0 = 0
p = 0.7
gamma = 0.3
maxT = 40
numSims = 100

#Set random seed for reproducability
np.random.seed(1)

#Initialise individuals
population = []
for i in range(0, N):
    population.append(Individual([], "Uninitialised"))
population = formNetwork(population, "random")

#Initialise lists to hold infectious numbers for all runs
listofIs = []

#Run simulation multiple times
for n in range(0, numSims):
    
    #Initialise infectious numbers for this run
    I = np.zeros(maxT + 1)
    I[0] = I0

    #Initialise individual states
    #Make the first S0 individuals susceptible
    for i in range(0, N-I0-R0):
        population[i].state = "S"
    #Make the next I0 individuals infectious
    for i in range(N-I0-R0, N-R0):
        population[i].state = "I"
    #Make the next R0 individuals recovered
    for i in range(N-R0, N):
        population[i].state = "R"
    
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
                    #If neighbour is infectious, current individual is infected with probability p
                    if (currentNeighbour.state == "I"):
                        rand = np.random.uniform(0, 1)
                        if (rand <= p):
                            infected = True
                            #If current individual is infected by one neighbour, there is no need to check if they are infected by others
                            break
                if infected:
                    nextStates[i] = "I"
                else:
                    nextStates[i] = "S"
                
            #Determine next state for infectious individuals
            elif (currentIndividual.state == "I"):
                #Current individual recovers with probability gamma
                rand = np.random.uniform(0, 1)
                if (rand <= gamma):
                    nextStates[i] = "R"
                else:
                    nextStates[i] = "I"

            #Determine next state for recovered individuals
            elif (population[i].state == "R"):
                #Recovered individuals always stay recovered
                nextStates[i] = "R"

        #Update states and count number of infectious individuals
        infectiousCount = 0
        for i in range(0, N):
            population[i].state = nextStates[i]
            if (population[i].state == "I"):
                infectiousCount += 1
        I[t] = infectiousCount

    #Add infectious numbers for this run to list
    listofIs.append(I)
    print("Completed simulation " + str(n))

#Plot infectious numbers over time for all runs
plt.rcParams.update({'font.size': 14})
plt.figure()
#plt.suptitle('Number of infectious individuals over time in SIR ? network model')

for n in range(0, numSims):
    plt.plot(np.linspace(0, maxT, maxT + 1), listofIs[n], 'r-')
#Calculate and plot mean infectious numbers
meanI = np.zeros(maxT + 1)
for t in range(0, maxT + 1):
    meanI[t] = np.mean([listofIs[n][t] for n in range(0, numSims)])
plt.plot(np.linspace(0, maxT, maxT + 1), meanI, 'k-', linewidth = 2.5)
plt.xlabel('Time')
plt.ylabel('I')

plt.show()
