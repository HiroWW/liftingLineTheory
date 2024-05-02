import numpy as np
import math
import matplotlib.pylab as plt

# Configuration parameters
N = 10  # This can now be any integer >= 1
S = 5.00  # wing area (m^2)
AR = 5.0  # Aspect ratio
taper = 1.0  # Taper ratio
alpha_twist = 0.0  # Twist angle (degrees)
i_w = 10.0  # Wing setting angle (degrees)
a_2d = 6.8754  # lift curve slope (1/rad)
alpha_0 = 0.0  # zero-lift angle of attack (degrees)

# Calculated values
b = math.sqrt(AR * S)  # Wing span (m)
MAC = S / b  # Mean Aerodynamic Chord (m)
Croot = (1.5 * (1 + taper) * MAC) / (1 + taper + taper ** 2)  # Root chord (m)
theta = np.linspace(math.pi / (2 * N), math.pi / 2, N, endpoint=True)
alpha = np.linspace(i_w + alpha_twist, i_w, N)
z = (b / 2) * np.cos(theta)
c = Croot * (1 - (1 - taper) * np.cos(theta))
mu = c * a_2d / (4 * b)
LHS = mu * (alpha - alpha_0) / 57.3

# Compute RHS matrix
RHS = []
for i in range(1, 2 * N, 2):
    RHS_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(theta)))
    RHS.append(RHS_iter)

# Ensure the RHS matrix has N columns
RHS = np.array(RHS).T
if RHS.shape[1] > N:
    RHS = RHS[:, :N]

# Solve the system using pseudoinverse for numerical stability
inv_RHS = np.linalg.pinv(RHS)
ans = np.matmul(inv_RHS, LHS)

# Calculate CL at each section
mynum = (4 * b) / c
CL = sum((np.sin((i+1) * theta)) * ans[i] * mynum for i in range(N))

# Append 0 at the beginning of CL array for plotting
CL1 = np.append(0, CL)
y_s = np.append(b / 2, z)

# Calculate and print the overall CL of the wing
CL_wing = math.pi * AR * ans[0]
print(f"Overall CL of the wing: {CL_wing:.2f}")

# Plotting
plt.plot(y_s, CL1, marker="o")
plt.title("Lifting Line Theory\nElliptical Lift Distribution")
plt.xlabel("Semi-span location (m)")
plt.ylabel("Lift coefficient")
plt.grid()
plt.show()


