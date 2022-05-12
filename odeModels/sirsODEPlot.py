import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Set parameters
N = 1000
beta = 0.5
gamma = 1
delta = 0.5
S0 = 500
I0 = 500
#Range of time values for which solution will be calculated
plotT = np.linspace(0, 40, 82)

#Set up system of equations and initial conditions (here, X = [S, I])
def dX_dt(X, t):
    return [delta*N - delta*X[0] - delta*X[1] - beta*X[1]*X[0]/N,
            beta*X[1]*X[0]/N - gamma*X[1]]

X0 = [S0, I0]

#Solve the system numerically
sol = odeint(dX_dt, X0, plotT)
S = sol[:,0]
I = sol[:,1]
R = N - S - I

#Create phase portrait
#Create grid
plotSArray = np.linspace(0, N, 16)
plotIArray = np.linspace(0, N, 16)
plotS, plotI = np.meshgrid(plotSArray, plotIArray)
#Initialise 2D arrays for horizontal and vertical components
u = np.zeros(plotS.shape)
v = np.zeros(plotI.shape)
#Populate arrays
for i in range(0, len(plotSArray)):
    for j in range(0, len(plotIArray)):
        #Get current coordinate
        X = [plotS[i, j], plotI[i, j]]
        #Calculate vector components (using t = 0 since system is autonomous)
        u[i, j] = dX_dt(X, 0)[0]
        v[i, j] = dX_dt(X, 0)[1]


#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()

#Plot phase portrait
plt.subplot(1,2,1)
plt.quiver(plotS, plotI, u, v)
plt.plot(S, I, 'k-')
plt.plot(np.linspace(0, N, 11), N- np.linspace(0, N, 11), 'r-')
plt.xlabel('S')
plt.ylabel('I')
#Plot numerical solution
plt.subplot(1,2,2)
plt.plot(plotT, I, 'r-', plotT, S, 'b-', plotT, R, 'g-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.ylim([0, N])
plt.legend(['Infectious', 'Susceptible', 'Recovered'])
plt.text(10, 550, r'$\beta = $' + str(beta) + ', $\gamma = $' + str(gamma) + ', $\delta = $' + str(delta) + '$, S_{0} = $' + str(S0) + ', $I_{0} = $' + str(I0))

plt.show()
