import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Set parameters for two different models and define the corresponding functions (where X[0] = S and X[1] = I)
N = 1000
#Plot 1
beta1 = 2
gamma1 = 1
def dX1_dt(X, t):
    return [-beta1*X[1]*X[0]/N, beta1*X[1]*X[0]/N - gamma1*X[1]]
#Plot 2
beta2 = 0.5
gamma2 = 1
def dX2_dt(X, t):
    return [-beta2*X[1]*X[0]/N, beta2*X[1]*X[0]/N - gamma2*X[1]]


#Solve the systems numerically for specified initial conditions
plotT = np.linspace(0, 40, 82)
#Plot 1
X01 = [900, 100]
sol1 = odeint(dX1_dt, X01, plotT)
S1 = sol1[:,0]
I1 = sol1[:,1]
#Plot 2
X02 = [500, 500]
sol2 = odeint(dX2_dt, X02, plotT)
S2 = sol2[:,0]
I2 = sol2[:,1]

#Create phase portraits for systems
#Create grid
plotSArray = np.linspace(0, N, 21)
plotIArray = np.linspace(0, N, 21)
plotS, plotI = np.meshgrid(plotSArray, plotIArray)
#Initialise arrays for horizontal and vertical components
#Plot 1
u1 = np.zeros(plotS.shape)
v1 = np.zeros(plotI.shape)
#Plot 2
u2 = np.zeros(plotS.shape)
v2 = np.zeros(plotI.shape)
#Populate arrays
for i in range(0, len(plotSArray)):
    for j in range(0, len(plotIArray)):
        #Get current coordinate
        X = [plotS[i, j], plotI[i, j]]
        #Retreive components (using t = 0 since system is autonomous)
        #Plot 1
        u1[i, j] = dX1_dt(X, 0)[0]
        v1[i, j] = dX1_dt(X, 0)[1]
        #Plot 2
        u2[i, j] = dX2_dt(X, 0)[0]
        v2[i, j] = dX2_dt(X, 0)[1]


#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()
plt.suptitle('Phase portraits and number of infectious, susceptible and recovered individuals over time in SIR ODE model')
#Plot phase portrait for system 1
plt.subplot(2,2,1)
plt.quiver(plotS, plotI, u1, v1)
plt.plot(S1, I1, 'k-')
plt.plot(np.linspace(1, 1000, 11), 1000- np.linspace(1, 1000, 11), 'r-')
plt.xlabel('S')
plt.ylabel('I')
#Plot numerical solution of system 1
plt.subplot(2,2,2)
plt.plot(plotT, I1, 'r-', plotT, S1, 'b-', plotT, N - S1 - I1, 'g-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible', 'Recovered'])
plt.text(20, 700, r'$\beta = $' + str(beta1) + ', $\gamma = $' + str(gamma1) + '$, S_{0} = $' + str(X01[0]) + ', $I_{0} = $' + str(X01[1]))

#Plot phase portrait for system 2
plt.subplot(2,2,3)
plt.quiver(plotS, plotI, u2, v2)
plt.plot(S2, I2, 'k-')
plt.plot(np.linspace(1, 1000, 11), 1000- np.linspace(1, 1000, 11), 'r-')
plt.xlabel('S')
plt.ylabel('I')
#Plot numerical solution of system 2
plt.subplot(2,2,4)
plt.plot(plotT, I2, 'r-', plotT, S2, 'b-', plotT, N - S2 - I2, 'g-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible', 'Recovered'], loc = 'lower left', bbox_to_anchor = (0.7, 0.1))
plt.text(20, 300, r'$\beta = $' + str(beta2) + ', $\gamma = $' + str(gamma2) + '$, S_{0} = $' + str(X02[0]) + ', $I_{0} = $' + str(X02[1]))

plt.show()
