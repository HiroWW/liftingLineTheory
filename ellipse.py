import numpy as np
import math
import matplotlib.pylab as plt

# Configuration parameters
N = 50  # 分割数
S = 5.00  # 翼面積 (m^2)
AR = 5.0  # アスペクト比
alpha_twist = 0.0  # ツイスト角 (degrees)
i_w = 10.0  # 翼の取り付け角 (degrees)
a_2d = 2 * math.pi  # 二次元の揚力曲線の傾き (1/rad)
alpha_0 = 0.0  # ゼロ揚力迎角 (degrees)

# Calculated values
b = math.sqrt(AR * S)  # 翼幅 (m)
theta = np.linspace(0.01, np.pi - 0.01, N)  # 翼のセグメントの中心角
z = (b / 2) * np.cos(theta)  # 翼の半径方向の位置
c = 4 * S / (np.pi * b) * np.sqrt(1 - (2 * z / b)**2)  # 楕円翼の翼弦長の分布
alpha = i_w + alpha_twist - (alpha_twist / (np.pi/2)) * theta  # 各セグメントの迎角
mu = c * a_2d / (4 * b)  # 揚力線の無次元パラメータ

# Compute LHS vector
LHS = mu * (alpha - alpha_0) * (np.pi/180)  # radianへ変換

# Compute RHS matrix
RHS = []
for i in range(1, 2 * N + 1, 2):
    sin_theta = np.sin(theta)
    sin_theta[sin_theta == 0] = 1e-10  # sin(theta)が0になるのを防ぐための小さな値を設定
    RHS_iter = np.sin(i * theta) * (1 + (mu * i) / sin_theta)
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
CL = np.array([sum(np.sin((2 * i + 1) * theta) * ans[i] * mynum[i]) for i in range(N)])

# Append 0 at the beginning and end of CL array for plotting
CL1 = np.append(0, np.append(CL, 0))
y_s = np.append(-b / 2, np.append(z, b / 2))

# Calculate and print the overall CL of the wing
CL_wing = math.pi * AR * ans[0]
print(f"Prandtl Average CL: {CL_wing:.4f}")
print("Helmbold Eq     CL:", a_2d / (1 + a_2d / (math.pi * AR)) * i_w / 180 * math.pi)

# 計算されたCLを左右の翼に拡張する
CL_full = np.concatenate((CL1[::-1], CL1))
# 翼の全スパンの位置を計算する
y_s_full = np.concatenate((-y_s[::-1], y_s))

# Plotting
plt.plot(y_s_full, CL_full, marker="o")
plt.title("Lifting Line Theory\nFull Span Lift Distribution")
plt.xlabel("Span location (m)")
plt.ylabel("Lift coefficient")
plt.grid()
plt.show()
