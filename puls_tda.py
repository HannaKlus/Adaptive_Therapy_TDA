
#-LIBRARIES-
import numpy as np
import matplotlib.pyplot as plt
from ripser import ripser
from persim import plot_diagrams

#-TIME AND SPACE CONFIGURATION-
dt = 0.01
t_max = 50
t_steps = np.arange(0, t_max, dt)
N = 100
L = 400
dx = L/N
x = np.linspace(0, L, N)
Nt = len(t_steps)
t_mutation = 35
med_effect = 1.1

#-INITIALIZATION-
v = np.zeros(N)
tumor_weight_history = np.zeros(Nt)

#-MODEL PARAMETERS-
#Diffusion constant
D_v = 0.1
#Growth rates
r_v = 1.0 #Tumor proliferation rate

#-INITIAL CONDITIONS (t=0)
middle = N // 2
width = 10
v[middle-width : middle+width] = 0.1

#-NUMERICAL INTEGRATION (Euler Method)-
for i in range(Nt):
    if t_steps[i] > t_mutation:
        med_effect = 0.1
    else:
        if t_steps[i] % 6 < 3:
            med_effect = 1.1
        else:
            med_effect = 0.0

    current_weight = np.sum(v) * dx
    tumor_weight_history[i] = current_weight

    #Tumor Diffusion
    lap_v = np.zeros(N)
    lap_v[1:-1] = (v[:-2] - 2 * v[1:-1] + v[2:]) / dx ** 2
    #
    dv = r_v * v * (1 - v) - med_effect * v
    #
    v_new = v + (dv + lap_v * D_v) * dt
    #
    v_new[v_new < 0.0] = 0.0
    #
    v = v_new

#TDA - Taken's theorem and Time delayed embedding
tau = 22
X_tda = []

for i in range(len(tumor_weight_history) - tau):
    x_now = tumor_weight_history[i]
    x_delayed = tumor_weight_history[i + tau]
    X_tda.append([x_now, x_delayed])
X_tda = np.array(X_tda)

window_size = 1000
step_size = 50
h1_persistence = []
time_axis = []

for start in range(0, len(tumor_weight_history) - tau - window_size, step_size):
    end = start + window_size
    window = X_tda[start:end]

    result = ripser(window, maxdim=1)
    diagrams = result['dgms']

    if len(diagrams[1]) == 0:
        h1_persistence.append(0.0)

    else:
        temorary_persistence = []
        for pt in diagrams[1]:
            persistence = pt[1] - pt[0]
            temorary_persistence.append(persistence)
        h1_persistence.append(np.max(temorary_persistence))
    time_axis.append(t_steps[start + window_size//2])

plt.figure(figsize=(10, 5))
plt.plot(time_axis,h1_persistence, color='red', linewidth=2)
plt.title('Mutation Detector')
plt.xlabel('Time [a.u.]')
plt.ylabel('Max H1 Persistence')
plt.grid(True, alpha=0.3)
plt.show()





