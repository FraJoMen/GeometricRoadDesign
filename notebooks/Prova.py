import numpy as np
import matplotlib.pyplot as plt
from scipy.special import fresnel

def clotoide(t, A):
    S, C = fresnel(t)
    x = A * C
    y = A * S
    return x, y

def versore_tangente(t, A):
    dx_dt = A * np.cos(np.pi * t**2 / 2)
    dy_dt = A * np.sin(np.pi * t**2 / 2)
    T = np.vstack((dx_dt, dy_dt))
    T_norm = T / np.linalg.norm(T, axis=0)
    return T_norm

def versore_normale(tangente):
    N = np.vstack((-tangente[1], tangente[0]))
    return N

# Parametri
A = 1
R = 0.8
L = A/(R*np.pi)

t = np.linspace(0, L/A, 100)

# Genera i punti della clotoide
x, y = clotoide(t, A)

# Calcola i versori tangenti e normali in punti specifici
t_end = L/A
T_end = versore_tangente(np.array([t_end]), A)
N_end = versore_normale(T_end)
x_end, y_end = clotoide(np.array([t_end]), A)

# Calcola il centro del cerchio osculatore
center_x = x_end + R * N_end[0]
center_y = y_end + R * N_end[1]

# Genera i punti del cerchio osculatore
theta = np.linspace(0, 2*np.pi, 100)
circle_x = center_x + R * np.cos(theta)
circle_y = center_y + R * np.sin(theta)

# Plot della clotoide e del cerchio osculatore
plt.plot(x, y, label='Clotoide')
plt.plot(circle_x, circle_y, 'g--', label='Cerchio Osculatore')
plt.plot(center_x, center_y, 'go', label='Centro del Cerchio Osculatore')
plt.quiver(x_end, y_end, T_end[0], T_end[1], color='r', scale=20, label='Versore Tangente')
plt.quiver(x_end, y_end, N_end[0], N_end[1], color='b', scale=20, label='Versore Normale')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Clotoide con Cerchio Osculatore')
#plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()
