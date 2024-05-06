import numpy as np
import math
import matplotlib.pylab as plt

# Configuration parameters
N = 50  # 分割数
S = 5.00  # 翼面積 (m^2)
AR = 5.0  # アスペクト比
taper = 1.0  # テーパー比
alpha_twist = 0.0  # ツイスト角 (degrees)
i_w = 10.0  # 翼の取り付け角 (degrees)
a_2d = 2 * math.pi  # 二次元の揚力曲線の傾き (1/rad)
alpha_0 = 0.0  # ゼロ揚力迎角 (degrees)

# Calculated values
b = math.sqrt(AR * S)  # 翼幅 (m)
MAC = S / b  # 平均空力翼弦 (m)
Croot = (1.5 * (1 + taper) * MAC) / (1 + taper + taper ** 2)  # 翼根のコード長 (m)
theta = np.linspace(math.pi / (2 * N), math.pi / 2, N, endpoint=True)  # 翼のセグメントの中心角
alpha = np.linspace(i_w + alpha_twist, i_w, N)  # 各セグメントの迎角
z = (b / 2) * np.cos(theta)  # 翼の半径方向の位置
c = Croot * (1 - (1 - taper) * np.cos(theta))  # 翼弦長の分布
mu = c * a_2d / (4 * b)  # 揚力線の無次元パラメータ

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
# plt.show()


import numpy as np
import matplotlib.pyplot as plt

def parse_ivs_result(filename):
    cl_values = []
    with open(filename, 'r') as file:
        for line in file:
            if "cl and cd is" in line:
                parts = line.strip().split()
                cl_value = parts[4]  # The cl value is the fifth word in the line
                cl_values.append(float(cl_value))
    return cl_values

def parse_divided_points(filename):
    x_values = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            x_values.append(float(parts[1]))  # Second column values
    return x_values

def generate_ellipse(x_range, AR, alpha):
    alpha_rad = alpha * np.pi / 180  # Convert alpha from degrees to radians
    a = 2 * np.pi / (1 + 2 * np.pi / (np.pi * AR))
    target_mean_y = a * alpha_rad  # Desired mean value of the y-values of the ellipse

    # b is initially 1, and we adjust it to scale the average to the target_mean_y
    y_values = np.sqrt(1 - (x_range / max(x_range))**2)  # Base ellipse with b = 1
    current_mean_y = np.mean(y_values)  # Calculate current mean of base ellipse y-values
    
    # Scale y_values so their mean matches target_mean_y
    scale_factor = target_mean_y / current_mean_y
    y_values *= scale_factor  # Scale the ellipse to match the desired mean

    return y_values

def plot_data(x_values, cl_values, x_range, ellipse_y_values):
    # plt.figure(figsize=(10, 6))
    plt.plot(x_values, cl_values, marker='o', linestyle='-', label='pointsource Cl')
    plt.plot(x_range, ellipse_y_values, linestyle='--', color='red', label='Ellipse curve')
    plt.xlabel('SPAN [m]')
    plt.ylabel('Cl')
    # plt.title('Plot of CL vs. Corresponding X Values and Ellipse Curve')
    plt.grid(True)
    plt.ylim(bottom=0)  # Start y-axis at 0
    plt.legend()
    plt.show()

# Load data
cl_values = parse_ivs_result('IVSRESULT.txt')
x_values = parse_divided_points('Divided_Points.txt')

# Calculate the x-axis range modifications
delta = (x_values[-1] - x_values[0]) / (len(x_values) - 1) / 2
x_values[0] -= delta
x_values[-1] += delta

# Generate a denser x range for ellipse
x_range = np.linspace(x_values[0], x_values[-1], num=1000)  # Create 1000 points for a smoother ellipse

# Generate ellipse values
ellipse_y_values = generate_ellipse(x_range, AR=5, alpha=10)

# Plotting the data
plot_data(x_values, cl_values, x_range, ellipse_y_values)

# Calculate and print the average of cl_values
average_cl = np.mean(cl_values)
modified_average_cl = np.mean(cl_values[1:-1]) * ( len(cl_values)-2 ) / len(cl_values)
print("Average of CL values:", average_cl)
print("modified:", modified_average_cl)
print("Theorotical CL:", np.mean(ellipse_y_values))