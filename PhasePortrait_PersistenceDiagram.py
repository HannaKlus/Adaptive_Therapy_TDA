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
R = 10 #initial tumor radius
V[dist < R**2] = 0.1
W[dist < R**2] = 0.005

#-NUMERICAL INTEGRATION (Euler Method)-
for i in range(Nt):
    if t[i] > t_mutation:
        med_effect = 0.1
    else:
        if t[i] % 6 < 3:
            med_effect = 1.1
        else:
            med_effect = 0.0
    # Current state
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
tau = 150
X_tda = []

for i in range(len(tumor_weight_history) - tau):
    x_now = tumor_weight_history[i]
    x_delayed = tumor_weight_history[i + tau]
    X_tda.append([x_now, x_delayed])
X_tda = np.array(X_tda)
#-PLOTTING-

plt.figure(figsize=(12, 5))


plt.subplot(1, 2, 1)
plt.plot(X_tda[:, 0], X_tda[:, 1], color='purple', alpha=0.6)
plt.title('1. Phase Portrait')
plt.xlabel('Mass (t)')
plt.ylabel(f'Mass (t + {tau})')
plt.grid(True, alpha=0.3)
plt.xlim(-10,650)
plt.ylim(-10,650)

plt.scatter(X_tda[0,0], X_tda[0,1], color='green', label='Start (Benign)')
plt.scatter(X_tda[-1,0], X_tda[-1,1], color='red', label='End (Invasive)')
plt.legend()

X_tda = X_tda[::20]
result = ripser(X_tda, maxdim=1)
diagrams = result['dgms']

plt.subplot(1, 2, 2)
plot_diagrams(diagrams, show=False)
plt.title('2. Topological print (Persistence Diagram)')
plt.tight_layout()
plt.savefig('Phase Portrait persistence diagram.png')
plt.show()




