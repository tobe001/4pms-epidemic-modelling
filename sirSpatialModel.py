import numpy as np
import matplotlib.pyplot as plt

#Create class for individuals to store their x and y coordinates and current state
class Individual:

    def __init__(self, x, y, state):
        self.x = x
        self.y = y
        self.state = state

#Create function to calculate distance between two individuals
def dist(i1, i2):
    return np.sqrt(pow(i1.x - i2.x, 2) + pow(i1.y - i2.y, 2))

#Define parameters
N = 1000
I0 = 1
R0 = 0
gamma = 0.3
maxT = 40
numSims = 100

#Define contact probability function
def p(d):
    k = 0.9
    m = 1
    return k*np.exp(-m * d)

#Set random seed for reproducability
np.random.seed(1)

#Initialise individuals
population = []
for i in range(0, N):
    #Draw x and y coordinates from normal distribution
    x = np.random.normal(0, 10)
    y = np.random.normal(0, 10)
    population.append(Individual(x, y, "Unitialised"))

#Initialise lists to hold susceptible, infectious and recovered numbers for all runs
listofSs = []
listofIs = []
listofRs = []

#Run simulation multiple times
for n in range(0, numSims):
    
    #Initialise susceptible, infectious and recovered numbers
    I = np.zeros(maxT + 1)
    S = np.zeros(maxT + 1)
    R = np.zeros(maxT + 1)
    I[0] = I0
    R[0] = R0
    S[0] = N - I0 - R0

    #Initialise the states of individuals
    #Make the first S0 susceptible
    for i in range(0, N-I0-R0):
        population[i].state = "S"
    #Make the next I0 infectious
    for i in range(N-I0-R0, N-R0):
        population[i].state = "I"
    #Make the nest R0 recovered
    for i in range(N-R0, N):
        population[i].state = "R"
    
    #Run simulation
    for t in range(1, maxT + 1):
        #Create list to hold states for next time step
        nextStates = [""] * N
        #For each individual:
        for i in range(0, N):
            currentIndividual = population[i]
            
            #Determine next state for susceptible individuals
            if (currentIndividual.state == "S"):
                infected = False
                #For each infectious individual:
                for j in range(0, len(population)):
                    currentNeighbour = population[j]
                    if (currentNeighbour.state == "I"):
                        #Infect with probability p(dij)
                        dij = dist(currentIndividual, currentNeighbour)
                        contactProb = p(dij)
                        rand = np.random.uniform(0, 1)
                        if (rand <= contactProb):
                            infected = True
                            #If current individual is infected once, no need to check if they are infected by any other individuals
                            break
                if infected:
                    nextStates[i] = "I"
                else:
                    nextStates[i] = "S"
                
            #Determine next state for infectious individuals
            elif (currentIndividual.state == "I"):
                #Recover with probability gamma
                rand = np.random.uniform(0, 1)
                if (rand <= gamma):
                    nextStates[i] = "R"
                else:
                    nextStates[i] = "I"

            #Determine next state for recovered individuals
            elif (population[i].state == "R"):
                #Always remain recovered
                nextStates[i] = "R"

        #Update states and count number of susceptible, infectious and recovered
        susceptibleCount = 0
        infectiousCount = 0
        recoveredCount = 0
        for i in range(0, N):
            population[i].state = nextStates[i]
            if (population[i].state == "S"):
                susceptibleCount += 1
            elif (population[i].state == "I"):
                infectiousCount += 1
            elif (population[i].state == "R"):
                recoveredCount += 1
        S[t] = susceptibleCount
        I[t] = infectiousCount
        R[t] = recoveredCount

    #Add susceptible, infectious and recovered numbers for this run to lists
    listofSs.append(S)
    listofIs.append(I)
    listofRs.append(R)
    print("Completed simulation " + str(n))

#Plot susceptible, infectious and recovered numbers over time for all runs
plt.rcParams.update({'font.size': 14})
plt.figure()
plt.suptitle('Number of susceptible, infectious and recovered individuals over time in SIR spatial model')
#Plot susceptible numbers
plt.subplot(1, 3, 1)
for n in range(0, numSims):
    plt.plot(np.linspace(0, maxT, maxT + 1), listofSs[n], 'c-')
#Calculate and plot mean susceptible numbers
meanS = np.zeros(maxT + 1)
for t in range(0, maxT + 1):
    meanS[t] = np.mean([listofSs[n][t] for n in range(0, numSims)])     
plt.plot(np.linspace(0, maxT, maxT + 1), meanS, 'k-', linewidth = 2.5)
plt.ylim([0, N])
plt.xlabel('Time')
plt.ylabel('S')

#Plot infectious numbers
plt.subplot(1, 3, 2)
for n in range(0, numSims):
    plt.plot(np.linspace(0, maxT, maxT + 1), listofIs[n], 'r-')
#Calculate and plot mean infectious numbers
meanI = np.zeros(maxT + 1)
for t in range(0, maxT + 1):
    meanI[t] = np.mean([listofIs[n][t] for n in range(0, numSims)]) 
plt.plot(np.linspace(0, maxT, maxT + 1), meanI, 'k-', linewidth = 2.5)
plt.ylim([0, N])
plt.xlabel('Time')
plt.ylabel('I')

#Plot recovered numbers
plt.subplot(1, 3, 3)
for n in range(0, numSims):
    plt.plot(np.linspace(0, maxT, maxT + 1), listofRs[n], 'g-')
#Calculate and plot mean recovered numbers
meanR = np.zeros(maxT + 1)
for t in range(0, maxT + 1):
    meanR[t] = np.mean([listofRs[n][t] for n in range(0, numSims)])    
plt.plot(np.linspace(0, maxT, maxT + 1), meanR, 'k-', linewidth = 2.5)
plt.ylim([0, N])
plt.xlabel('Time')
plt.ylabel('R')

plt.show()
