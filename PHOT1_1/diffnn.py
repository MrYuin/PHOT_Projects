import numpy as np


def diffnn(x):
    """
    计算数值微分

    参数:
        x: 输入数组

    返回:
        f: 微分结果数组
    """
    n = len(x)
    f = np.zeros(n)

    # 中间点使用中心差分
    f[1:n - 1] = 0.5 * (x[2:n] - x[0:n - 2])

    # 边界点使用前向/后向差分
    f[0] = x[1] - x[0]
    f[n - 1] = x[n - 1] - x[n - 2]

    return f