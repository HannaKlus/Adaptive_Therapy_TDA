
#-LIBRARIES-
import numpy as np
import matplotlib.pyplot as plt
from ripser import ripser
from persim import plot_diagrams

#-TIME AND SPACE CONFIGURATION-
dt = 0.002
t_max = 60
t = np.arange(0, t_max, dt)
Nx = 201
Ny = 201
L = 400
dx = L/(Nx-1)
dy = L/(Ny-1)
x = np.linspace(0, L, Nx)
y = np.linspace(0, L, Ny)
Nt = len(t)

#-INITIALIZATION-
U = np.zeros((Ny, Nx))
V = np.zeros((Ny, Nx))
W = np.zeros((Ny, Nx))
tumor_weight_history = np.zeros(Nt)

#-MODEL PARAMETERS-
#Diffusion constants
D_u = 0
D_v = 0.00004
D_w = 1.0
#Other parameters
t_mutation = 35
med_effect = 0.0
delta_1 = 12.5
delta_3 = 70.0
#Growth rates
r_v = 1.0 #Tumor proliferation rate

#-INITIAL CONDITIONS (t=0)
U[:, :] = 1.0
#
X, Y = np.meshgrid(x, y)
dist = (X - L/2)**2 + (Y - L/2)**2 #Distence from the center (matrix Nx x Ny)
R = 40 #initial tumor radius
V[dist < R**2] = 0.5
W[dist < R**2] = 0.005

M0 = np.sum(V) * dx**2 #Initial Mass
is_treated = True

#-NUMERICAL INTEGRATION (Euler Method)-
for i in range(Nt):
    current_mass = np.sum(V) * dx ** 2
    if t[i] > t_mutation:
        med_effect = 0.1
    else:
        if is_treated and current_mass <= 0.5 * M0:
            is_treated = False
        elif not is_treated and current_mass >= M0:
            is_treated = True
        med_effect = 1.1 if is_treated else 0.0


    #Acid Diffusion
    lap_w = np.zeros((Ny, Nx))
    lap_w[1:-1, 1:-1] = (W[:-2, 1:-1] + W[2:, 1:-1] + W[1:-1, 2:] + W[1:-1, :-2] - 4 * W[1:-1, 1:-1]) / dx**2
    #Tumor Diffusion
    D = D_v * (1 - U)

    dV_right = (V[1:-1, 2:] -V[1:-1, 1:-1])
    dV_left =  (V[1:-1, :-2] - V[1:-1, 1:-1])
    dV_up = (V[:-2, 1:-1] - V[1:-1, 1:-1])
    dV_down = (V[2:, 1:-1] - V[1:-1, 1:-1])

    D_right = (D[1:-1, 2:] + D[1:-1, 1:-1]) / 2.0
    D_left = (D[1:-1, :-2] + D[1:-1, 1:-1]) /2.0
    D_up = (D[:-2, 1:-1] + D[1:-1, 1:-1]) /2.0
    D_down = (D[2:, 1:-1] + D[1:-1, 1:-1]) /2.0

    lap_v = np.zeros((Ny, Nx))
    #Sum of fluxes
    lap_v[1:-1, 1:-1] = (D_right * dV_right + D_left * dV_left + D_up * dV_up + D_down  * dV_down) / dx**2

    du = U * (1 - U) - delta_1 * U * W
    dv = r_v * V * (1 - V) - med_effect * V
    dw = delta_3 * (V - W)
    #
    U = U + du * dt
    V = V + (dv + lap_v) * dt
    W = W + (dw + lap_w * D_w) * dt
    #
    U = np.clip(U, 0, 1)
    V = np.clip(V, 0, 1)
    W = np.clip(W, 0, None)

    #Neumann Boundary Conditions
    for matrix in (U, V, W):
        matrix[0, :] = matrix[1, :]
        matrix[-1, :] = matrix[-2, :]
        matrix[:, 0] = matrix[:, 1]
        matrix[:, -1] = matrix[:, -2]
    #
    tumor_weight_history[i] = np.sum(V) * dx**2


#TDA - Taken's theorem and Time delayed embedding
tau = 250
X_tda = []

for i in range(len(tumor_weight_history) - tau):
    x_now = tumor_weight_history[i]
    x_delayed = tumor_weight_history[i + tau]
    X_tda.append([x_now, x_delayed])
X_tda = np.array(X_tda)

window_size = 4000
step_size = 100
h1_persistence = []
time_axis = []

for start in range(0, len(tumor_weight_history) - tau - window_size, step_size):
    end = start + window_size
    window = X_tda[start:end]

    window = window[::20]

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
    time_axis.append(t[start + window_size//2])

plt.figure(figsize=(18, 6))

plt.subplot(1,2,1)
plt.plot(time_axis,h1_persistence, color='red', linewidth=2)
plt.title('Mutation Detector')
plt.xlabel('Time [a.u.]')
plt.ylabel('Max H1 Persistence')
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
plt.plot(t, tumor_weight_history, color='purple', label='Tumor Weight')
plt.tight_layout()
plt.title('Tumor Mass Dynamics')
plt.xlabel('Time [a.u.]')
plt.ylabel('Mass [a.u.]')
plt.grid(True, alpha=0.3)
plt.subplots_adjust(wspace=0.2)
plt.show()





