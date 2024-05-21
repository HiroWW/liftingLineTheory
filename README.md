# Lifting Line Theory
揚力線理論による主翼解析

## 実行方法

1. (1回目のみ) venv環境の作成 
    - `python3 -m venv .venv`
    - `source .venv/bin/activate`(bash)
    - `pip install -r requirements.txt`
2. (2回目以降) python仮想環境のアクティベート
    - `./.venv/Scripts/activate`
3. 実行
    - `python main.py`
4. 仮想環境の終了
    - `deactivate`

## 設定

翼の平面系をスクリプトの冒頭で記述する。各変数の意味は以下の通り。
``` python
N: int = 9  # 翼を分割するセグメントの数-1
S: float = 5.00  # 翼面積 (平方メートル)
AR: float = 5.0  # アスペクト比（翼幅の二乗を翼面積で割った値）
taper: float = 1.0  # テーパー比（翼端のコード長を翼根のコード長で割った値）
alpha_twist: float = 0.0  # ツイスト角（度）
i_w: float = 10.0  # 翼の取り付け角（度）
a_2d: float = 6.8754  # 二次元揚力線の傾き（1/ラジアン）
alpha_0: float = 0.0  # ゼロ揚力迎角（度）
```

## 解説

揚力線理論（Lifting Line Theory）は、主翼の揚力分布を解析するための方法です。この理論は、主翼を一次元の揚力線と見做して、その周りの流れを解析します。ここでは、揚力線理論に基づいたコードの解説を数式を交えて行います。

### 理論の概要

揚力線理論では、主翼を翼幅方向に長い直線（揚力線）としてモデル化し、この線から発生する渦（渦糸）が周囲の流れにどのように影響を与えるかを計算します。この理論は主に以下のステップで進行します：

1. **翼のジオメトリ設定**：翼の形状、アスペクト比（AR）、翼面積（S）、テーパー比（taper）などを設定します。
2. **離散化**：翼を複数のセグメントに分割し、それぞれのセグメントの中心での流れの条件を計算します。
3. **フーリエ級数展開**：揚力分布をフーリエ級数（正弦波の和）で表現し、各波の振幅を求めることで全体の揚力分布を求めます。
4. **連立方程式の解法**：フーリエ級数の係数を求めるために、線形代数の手法を用いて連立方程式を解きます。

### コードの解説

#### 1. パラメータの設定

```python
N = 10  # 分割数
S = 5.00  # 翼面積 (m^2)
AR = 5.0  # アスペクト比
taper = 1.0  # テーパー比
alpha_twist = 0.0  # ツイスト角 (degrees)
i_w = 10.0  # 翼の取り付け角 (degrees)
a_2d = 6.8754  # 二次元の揚力曲線の傾き (1/rad)
alpha_0 = 0.0  # ゼロ揚力迎角 (degrees)
```

ここで設定するパラメータは、翼の物理的な特性と飛行条件を定義します。

#### 2. 基本的な翼の特性の計算

```python
b = math.sqrt(AR * S)  # 翼幅 (m)
MAC = S / b  # 平均空力翼弦 (m)
Croot = (1.5 * (1 + taper) * MAC) / (1 + taper + taper ** 2)  # 翼根のコード長 (m)
theta = np.linspace(math.pi / (2 * N), math.pi / 2, N, endpoint=True)  # 翼のセグメントの中心角
alpha = np.linspace(i_w + alpha_twist, i_w, N)  # 各セグメントの迎角
z = (b / 2) * np.cos(theta)  # 翼の半径方向の位置
c = Croot * (1 - (1 - taper) * np.cos(theta))  # 翼弦長の分布
mu = c * a_2d / (4 * b)  # 揚力線の無次元パラメータ
```

$ b $、$ MAC $、$ Croot $ は、翼の物理的寸法を計算します。$ \theta $、$ \alpha $、$ z $、$ c $ は、翼のセグメントごとの角度、迎角、位置、弦長を求めるために使用します。

#### 3. 方程式の構築と解法

```python
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
```

ここでは、$ \text{RHS} $ 行列を計算し、フーリエ級数の各項に対応する $ a_n $ 係数を求めます。これにより、翼全体の揚力分布を解析できます。疑似逆行列を使用して数値的安定性を向上させています。

#### 4. 揚力係数の計算とグラフ表示

```python
mynum = (4 * b) / c
CL = sum((np.sin((i+1) * theta)) * ans[i] * mynum for i in range(N))

# Append 0 at the beginning of CL array for plotting
CL1 = np.append(0, CL)
y_s = np.append(b / 2, z)

# Plotting
plt.plot(y_s, CL1, marker="o")
plt.title("Lifting Line Theory\nElliptical Lift Distribution")
plt.xlabel("Semi-span location (m)")
plt.ylabel("Lift coefficient")
plt.grid()
plt.show()
```

最後に、計算された係数を使用して各セグメントにおける揚力係数 $ C_L $ を計算し、翼の揚力分布をプロットします。