import numpy as np
import matplotlib.pyplot as plt

#Set parameters
beta = 1
gamma = 1.5
I0 = 400
N = 1000

#Range of time
plotT = np.linspace(0, 40, 82)

#Compute I(t) and S(t) for given t using the analytically solved equation
I = (beta - gamma)*N*I0/(((beta - gamma)*N - beta*I0)*np.exp(-(beta - gamma)*plotT) + beta*I0)
S = N - I

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()

plt.plot(plotT, I, 'r-', plotT, S, 'b-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible'])
plt.text(22, 750, r'$\beta = $' + str(beta) + ', $\gamma = $' + str(gamma) +  ', $I_{0} = $' + str(I0))

plt.show()
