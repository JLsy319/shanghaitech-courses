# Computational Biology Homework 9

## Q1

GROMACS力场

- 总能量函数:

$$
V(\mathbf{r}^N) = \sum_{\text{bonds}} \frac{1}{2} K_b (b - b_0)^2 + \sum_{\text{bonds}} \frac{1}{4} K_b' (b - b_0)^4 + \sum_{\text{angles}} \frac{1}{2} K_\theta (\cos\theta - \cos\theta_0)^2 \\

+ \sum_{\text{impropers}} \frac{1}{2} K_\xi (\xi - \xi_0)^2 + \sum_{\text{dihedrals}} K_\phi [1 + \cos(n\phi - \delta)] + \sum_{i<j}^{\text{non-bonded}} \left[ \frac{C_{12}(ij)}{r_{ij}^{12}} - \frac{C_6(ij)}{r_{ij}^6} + \frac{q_i q_j}{4\pi\epsilon_0\epsilon_r r_{ij}} \right]
$$

- 键伸缩：

$$
E_{\text{bond}} = \sum_{\text{bonds}} \frac{1}{2} K_b (b - b_0)^2 + \sum_{\text{bonds}} \frac{1}{4} K_b' (b - b_0)^4
$$

- 键角:

$$
E_{\text{angle}} = \sum_{\text{angles}} \frac{1}{2} K_\theta (\cos\theta - \cos\theta_0)^2
$$

- 二面角:

$$
E_{\text{improper}} = \sum_{\text{impropers}} \frac{1}{2} K_\xi (\xi - \xi_0)^2
$$

- 范德华相互作用:

$$
E_{\text{LJ}} = \sum_{i<j}^{\text{non-bonded}} \left[ \frac{C_{12}(ij)}{r_{ij}^{12}} - \frac{C_6(ij)}{r_{ij}^6} \right]
$$

- 静电相互作用:

$$
E_{\text{Coulomb}} = \sum_{i<j}^{\text{non-bonded}} \frac{q_i q_j}{4\pi\epsilon_0\epsilon_r r_{ij}}
$$

- Combining Rules:

$$
C_6(ij) = \sqrt{C_6(ii) \cdot C_6(jj)} \\
C_{12}(ij) = \sqrt{C_{12}(ii) \cdot C_{12}(jj)}
$$

## Q2

QSAR 方法最早由 Corwin Hansch 教授提出

他的主要贡献有：
1. 开创了 QSAR 领域的研究
2. 重新定义了分配系数，计算过程中引入了疏水性参数 $\pi$ 来解释药物的疏水效应
3. 运用多参数的回归分析，同时分析电子效应，立体效应和疏水效应对于生物活性的影响
4. 建立了包括数万种化合物的 QSAR 数据库
5. 将分子动力学模拟与 QSAR 相结合，在二氢叶酸还原酶方面的研究展示了 QSAR 和分子图形学在药物构效关系上的解释和应用

## Q3

首先选择一系列结构相似的化合物，然后找到他们的生物活性相关数据（单位通常是浓度），用 $\lg (1/C)$ 表示生物活性，此处 $C$ 的单位为摩尔浓度。  
然后计算 Molecular descriptor，通常包括疏水性参数 $\pi$，电子参数（例如 Hammett参数 $\sigma$）和立体参数（例如 Taft 立体参数 $E$ 或者摩尔折射率等）。  
之后建立数学模型，运用多参数回归分析的方法，建立生物活性与 Molecular descriptor 之间的关系，通常的线性 QSAR 关系如下：

$$
\lg (\frac{1}{C}) = a \cdot \text{疏水性参数} + b \cdot \text{电子参数} + c \cdot \text{立体参数} + k \text{（常数）}
$$

最后利用各种方式验证建立好的 QSAR 模型，验证完成后则可以用来预测新化合物的生物活性并指导药物设计


## Q4

SBDD: 需要拥有靶点的准确三维结构信息，尤其是清晰的活性位点构象和关键氨基酸残基的空间位置

LBDD: 需要拥有一系列已知生物活性的配体，且这一系列配体中要具有结构多样性和较大的生物活性跨度

## Q5

首先获取靶点的高分辨率三维结构，可以通过 X-ray, NMR 或者 Cryo-EM 等方法得到，或者可以通过同源建模这样的计算方法，根据结构详细分析靶点潜在的药物结合位点，包括关键的氨基酸残基，氢键和疏水区域等。  
然后通过 virtual screening 或者 de novo design 的方法，前者将化合物数据库中的分子与靶点进行分子对接，评估靶点和配体之间的亲和力；后者根据靶点活性位点的特征，从头进行设计。  
通过计算方法筛选出潜在的化合物后，利用靶点的三维结构信息进行优化，主要由分子对接，分子动力学模拟和分子图形学可视化分析这几种方法。  
最后对新得到的化合物进行生物活性测试，然后可以得到化合物与靶点结晶后的结构，进一步迭代化合物分子。

