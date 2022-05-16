import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Set parameters
N = 1000
I0 = 1
gamma = -np.log(1 - 0.3)
mu = -np.log(1 - 0.1)
muHat = -np.log(1 - 0.05)
n = 7
etaA = 2/3
etaB = 2/3
tauA = 2/3
tauB = 2/3
treatmentsAvailable = 500

maxT = 80
gammaRiskHats = np.linspace(0.3, 1.0, 71)

#Define functions determining how many treatments should be given to each group based on the total number of currently available treatments and the number of infectious individuals in each group
def fA(currentlyAvailable, IA, IB):
    return (IA / (IA + IB)) * currentlyAvailable

def fB(currentlyAvailable, IA, IB):
    return (IB / (IA + IB)) * currentlyAvailable

#Define functions to calculate number of treatments provided to each group from the number intended for each one as a function of the total number available and the number of infectious individuals in each group
def numTreatmentsAA(currentlyAvailable, IA, IB):
    return np.min([fA(currentlyAvailable, IA, IB), IA])

def numTreatmentsBB(currentlyAvailable, IA, IB):
    return np.min([fB(currentlyAvailable, IA, IB), IB])

def numTreatmentsBA(currentlyAvailable, IA, IB):
    return np.min([IA - numTreatmentsAA(currentlyAvailable, IA, IB), fB(currentlyAvailable, IA, IB) - numTreatmentsBB(currentlyAvailable, IA, IB)])

def numTreatmentsAB(currentlyAvailable, IA, IB):
    return np.min([IB - numTreatmentsBB(currentlyAvailable, IA, IB), fA(currentlyAvailable, IA, IB) - numTreatmentsAA(currentlyAvailable, IA, IB)])

#Calculate beta parameters
betaAA = etaA * tauA * n
betaAB = etaA * tauB * n
betaBA = etaB * tauA * n
betaBB = etaB * tauB * n

#Create list for storing final fatality numbers for each number of available treatments
finalDs = []

for gammaRiskHat in gammaRiskHats:

    gammaHat = -np.log(1 - gammaRiskHat)
    
    #Set up equations (where X = [SA, SB, IA, IB, TA, TB, R, D])
    def dX_dt(X, t):
        return [-betaAA*X[2]*X[0]/N - betaBA*X[3]*X[0]/N,
                -betaAB*X[2]*X[1]/N - betaBB*X[3]*X[1]/N,
                betaAA*X[2]*X[0]/N + betaBA*X[3]*X[0]/N - numTreatmentsAA(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) - numTreatmentsBA(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) - gamma*X[2] - mu*X[2],
                betaAB*X[2]*X[1]/N + betaBB*X[3]*X[1]/N - numTreatmentsBB(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) - numTreatmentsAB(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) - gamma*X[3] - mu*X[3],
                numTreatmentsAA(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) + numTreatmentsBA(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) - gammaHat*X[4] - muHat*X[4],
                numTreatmentsBB(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) + numTreatmentsAB(treatmentsAvailable - (X[4] + X[5]), X[2], X[3]) - gammaHat*X[5] - muHat*X[5],
                gamma*X[2] + gamma*X[3] + gammaHat*X[4] + gammaHat*X[5],
                mu*X[2] + mu*X[3] + muHat*X[4] + muHat*X[5]]

    #Specify initial conditions
    SA0 = (N - I0)/2
    SB0 = (N - I0)/2
    IA0 = I0/2
    IB0 = I0/2
    TA0 = 0
    TB0 = 0
    R0 = 0
    D0 = 0
    X0 = [SA0, SB0, IA0, IB0, TA0, TB0, R0, D0]
    
    #Solve equations over specified time scale and with specified initial conditions
    plotT = range(0, maxT + 1)
    sol = odeint(dX_dt, X0, plotT)
    S = sol[:,0] + sol[:,1]
    I = sol[:,2] + sol[:,3]
    T = sol[:,4] + sol[:,5]
    R = sol[:,6]
    D = sol[:,7]

    #Check that final number of deceased individuals appears to be converging. If it is not, we should simulate for longer.
    if (abs(D[maxT] - D[maxT - 1]) > 0.1):
        print("WARNING: Possible non-convergence of fatality numbers when T-Hat = " + str(treatmentsAvailable) + ". Consider increasing maxT.")

    #Add final fatlity numbers to list
    finalDs.append(D[maxT])
    print("Completed simulation for gammaRiskHat = " + str(gammaRiskHat))

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()
#Plot final fatality numbers against v
plt.plot(gammaRiskHats, finalDs, 'k-', linewidth = 2.5)
plt.xlabel(r'$\hat{\gamma_{risk}}$')
plt.ylabel('Final number of deceased individuals')
plt.xlim([0.3,1.0])
plt.ylim([0, 300])

plt.show()
