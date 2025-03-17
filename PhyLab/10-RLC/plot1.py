import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.optimize import fsolve

# 假设我们有 x 和 y 数据
x = np.array([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300])
y = np.array([32.4,64,104,156,232,380,606,696,552,408,320,256,224,208,184,168,132,124,120,108,102,88,82,74,72])

# 创建光滑曲线
x_smooth = np.linspace(x.min(), x.max(), 500)
spl = make_interp_spline(x, y, k=3)
y_smooth = spl(x_smooth)

# 找到最大值及其位置
max_y = np.max(y_smooth)
max_x = x_smooth[np.argmax(y_smooth)]

# 计算 y = max_y / sqrt(2) 的水平线
level = max_y / np.sqrt(2)

# 定义函数来找到交点
def find_intersections(x_smooth, y_smooth, level):
    intersections = []
    for i in range(len(x_smooth) - 1):
        if (y_smooth[i] - level) * (y_smooth[i + 1] - level) < 0:
            x_intersect = fsolve(lambda x: spl(x) - level, x_smooth[i])
            intersections.append(x_intersect[0])
    return intersections

# 找到交点
intersections = find_intersections(x_smooth, y_smooth, level)

# 绘制曲线和交点
plt.plot(x_smooth, y_smooth, label='Smooth Curve')
plt.plot(max_x, max_y, 'ro', label='Maximum Point')
plt.axhline(level, color='r', linestyle='--', label=f'y = {level:.2f}')
for intersect in intersections:
    plt.plot(intersect, level, 'go')

plt.legend()
plt.show()

print("Maximum point: (", max_x, ",", max_y, ")")
print("Intersections with level line:", intersections)

