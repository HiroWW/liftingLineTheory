import numpy as np
import math
import matplotlib.pylab as plt
import os
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Lifting Line Theory')
parser.add_argument('--out', type=str, required=True, help='output directory')
# optional argument for switching to different mode
parser.add_argument('--taper', action='store_true', help='use it if you want to calculate tapered wing.')
parser.add_argument('--ellipse', action='store_true', help='use it if you want to calculate ellipse wing.')
args = parser.parse_args()

# Configuration parameters
N = 50  # 分割数
S = 5.00  # 翼面積 (m^2)
AR = 5.0  # アスペクト比
alpha_twist = 0.0  # ツイスト角 (degrees)
i_w = 10.0  # 翼の取り付け角 (degrees)
a_2d = 2 * math.pi  # 二次元の揚力曲線の傾き (1/rad)
alpha_0 = 0.0  # ゼロ揚力迎角 (degrees)

# Configuration parameters for tapered wing
if args.taper:
    taper = 1.0  # テーパー比
    b = math.sqrt(AR * S)  # 翼幅 (m)
    MAC = S / b  # 平均空力翼弦 (m)
    Croot = (1.5 * (1 + taper) * MAC) / (1 + taper + taper ** 2)  # 翼根のコード長 (m)
    theta = np.linspace(math.pi / (2 * N), math.pi / 2, N, endpoint=True)  # 翼のセグメントの中心角
    alpha = np.linspace(i_w + alpha_twist, i_w, N)  # 各セグメントの迎角
    z = (b / 2) * np.cos(theta)  # 翼のスパン方向位置
    c = Croot * (1 - (1 - taper) * np.cos(theta))  # 翼弦長の分布
    mu = c * a_2d / (4 * b)  # 揚力線の無次元パラメータ
# Configuration parameters for ellipse wing
if args.ellipse:
    b = math.sqrt(AR * S)  # 翼幅 (m)
    MAC = S / b  # 平均空力翼弦 (m)
    Croot = 2  # 翼根のコード長 (m)
    theta = np.linspace(math.pi / (2 * N), math.pi / 2, N, endpoint=True)  # 翼のセグメントの中心角
    alpha = np.linspace(i_w + alpha_twist, i_w, N)  # 各セグメントの迎角
    z = (b / 2) * np.cos(theta)  # 翼のスパン方向位置
    c = Croot * np.sin(theta)  # 翼弦長の分布
    mu = c * a_2d / (4 * b)  # 揚力線の無次元パラメータ

# other way to define the wing
# c = f(z)
# just define the function f(z) and use it to calculate c

# plot c - z and show it by useing matplotlig
plt.plot(z, c)
plt.title("Chord Distribution")
plt.xlabel("Span location (m)")
plt.ylabel("Chord length (m)")
plt.grid()
plt.show()

# Compute LHS vector
LHS = mu * (alpha - alpha_0) / 57.3

# Compute RHS matrix
RHS = []
for i in range(1, 2 * N + 1, 2):
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
CL = sum((np.sin(( 2 * i + 1) * theta)) * ans[i] * mynum for i in range(N))

# Append 0 at the beginning of CL array for plotting
CL1 = np.append(0, CL)
y_s = np.append(b / 2, z)

# Calculate and print the overall CL of the wing
CL_wing = math.pi * AR * ans[0]
print(f"Prandtl Average CL: {CL_wing:.4f}")
print( "Helmbold Eq     CL:", a_2d / (1+a_2d/(math.pi*AR)) * i_w/180*math.pi)

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