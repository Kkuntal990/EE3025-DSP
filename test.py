#%%
import matplotlib.pyplot as plt
import numpy as np
import math
#%%
N = 8
xn = np.array([1, 2, 3, 4, 2, 1,0,0])
hn = np.zeros(N)
for i in range(N):
    hn[i] += (-0.5)**i
    if i >= 2:
        hn[i] += (-0.5)**(i - 2)


def rfft(a):
    n = a.size
    if n == 1:
        return a
    i = 1j
    w_n = np.exp(-2 * i * np.pi / float(n))
    w = 1
    a_0 = np.zeros(int(math.ceil(n / 2.0)), dtype=np.complex_)
    a_1 = np.zeros(n - a_0.shape[0], dtype=np.complex_)
    for index in range(0, n):
        if index % 2 == 0:
            a_0[index // 2] = a[index]
        else:
            a_1[index // 2] = a[index]
    y_0 = rfft(a_0)
    y_1 = rfft(a_1)
    y = np.zeros(n, dtype=np.complex_)
    for k in range(0, n // 2):
        y[k] = y_0[k] + w * y_1[k]
        y[k + n // 2] = y_0[k] - w * y_1[k]
        w = w * w_n
    return y

#%%
X = rfft(xn)
plt.title("x(n)")
plt.stem(xn)
plt.show()

plt.stem(np.abs(X))
plt.show()
plt.stem(np.angle(X, True))
plt.show()

# %%
