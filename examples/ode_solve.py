import numpy as np
import matplotlib.pyplot as plt

# Lorenz system derivatives
def lorenz_system(current_state, t):
    x, y, z = current_state
    sigma = 100.0
    rho = 2800.0
    beta = 4/3
    
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return np.array([dxdt, dydt, dzdt])

# Runge-Kutta 4th Order Integrator
def rk4_step(func, y, t, dt):
    k1 = func(y, t)
    k2 = func(y + 0.5 * dt * k1, t + 0.5 * dt)
    k3 = func(y + 0.5 * dt * k2, t + 0.5 * dt)
    k4 = func(y + dt * k3, t + dt)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

# Simulation parameters
dt = 0.001
num_steps = 100000000
initial_state = np.array([10.0, 1.0, 105.0])

# Initialize arrays to store results
states = np.zeros((num_steps + 1, 3))
states[0] = initial_state

# Integration loop
for i in range(num_steps):
    states[i+1] = rk4_step(lorenz_system, states[i], i * dt, dt)

# Plotting
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Extract x, y, z coordinates
x, y, z = states[:, 0], states[:, 1], states[:, 2]

# Plot the trajectory with a color gradient for time
ax.plot(x, y, z, lw=0.5, color='royalblue')
ax.set_title("Lorenz Attractor (RK4 Solver)")
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")

plt.show()