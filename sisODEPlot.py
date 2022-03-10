import numpy as np
import matplotlib.pyplot as plt

#Set parameters for two different plots
#Plot1
beta1 = 1.5
gamma1 = 1
I01 = 1
#Plot 2
beta2 = 1
gamma2 = 1.5
I02 = 400

N = 1000

#Compute I(t) and S(t) for t from 0 to 40 for both plots using the analytically solved equation
plotT = np.linspace(0, 40, 82)
#Plot 1
I1 = (beta1 - gamma1)*N*I01/(((beta1 - gamma1)*N - beta1*I01)*np.exp(-(beta1 - gamma1)*plotT) + beta1*I01)
S1 = N - I1
#Plot 2
I2 = (beta2 - gamma2)*N*I02/(((beta2 - gamma2)*N - beta2*I02)*np.exp(-(beta2 - gamma2)*plotT) + beta2*I02)
S2 = N - I2

#Plot results
plt.rcParams.update({'font.size': 14})
plt.figure()
plt.suptitle('Number of infectious and susceptible individuals over time in SIS ODE model')
#Plot 1
plt.subplot(1, 2, 1)
plt.plot(plotT, I1, 'r-', plotT, S1, 'b-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible'])
plt.text(22, 750, r'$\beta = $' +str(beta1) + ', $\gamma = $' + str(gamma1) +  ', $I_{0} = $' + str(I01))
#Plot 2
plt.subplot(1, 2, 2)
plt.plot(plotT, I2, 'r-', plotT, S2, 'b-')
plt.xlabel('Time')
plt.ylabel('Number of individuals')
plt.legend(['Infectious', 'Susceptible'])
plt.text(22, 700, r'$\beta = $' + str(beta2) + ', $\gamma = $' + str(gamma2) + ', $I_{0} = $' + str(I02))

plt.show()
