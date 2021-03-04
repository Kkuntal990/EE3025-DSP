#%%
import matplotlib.pyplot as plt
import numpy as np

import timeit

#If using termux
import subprocess
import shlex
from numpy.core.fromnumeric import mean, size
from numpy.lib.function_base import average
#end if

#%%
N = 8
xn = np.array([1,2,3,4,2,1,0,0])
hn = np.zeros(N)
for i in range(N):
    hn[i] += (-0.5)**i
    if i >= 2:
        hn[i] += (-0.5)**(i - 2)
#%%
def DTFT(x, inverse=-1):
    M = x.shape[0]
    if M & (M-1) != 0:
        tmp = 1 << (M-1).bit_length()
        x = np.pad(x, [0, tmp - M], mode='constant')
        M = x.shape[0]

    W = np.zeros((M, M), dtype=complex)

    for n in range(M): 
        for k in range(M):
            W[n][k] = np.exp(inverse * 2j * np.pi * k * n / M)

    return np.matmul(W, x)

def FFT(x,inverse=-1):
    n = x.shape[0]
    if n == 1:
        return x
    if n&(n-1) != 0:
        tmp = 1 << (n-1).bit_length()
        x = np.pad(x, [0,tmp - n], mode='constant')
        n = x.shape[0]
    X1 = FFT(x[::2], inverse)
    X2 = FFT(x[1::2], inverse)
    w = np.exp(inverse*2j*np.pi*np.arange(n)/n)
    
    return np.concatenate([X1 + w[:n//2]*X2, X1 + w[n//2:]*X2])

def IFFT(X):
    n = X.shape[0]
    if n & (n-1) != 0:
        tmp = 1 << (n-1).bit_length()
        X = np.pad(x, [0, tmp - n], mode='constant')
        n = X.shape[0]
    return (1/n)*(FFT(X,inverse = 1))

#%%
X = FFT(xn)

plt.figure(0)
plt.stem(xn)
plt.xlabel('$n$')
plt.ylabel('$x(n)$')
plt.grid()
plt.show()
plt.savefig('../figs/x_n.pdf')
plt.savefig('../figs/x_n.eps')
#subprocess.run(shlex.split("termux-open ../figs/x_n.pdf"))  #if using termex

plt.figure(1,figsize=(10,4))

plt.subplot(1,2,1)
plt.stem(np.abs(X))
plt.grid()
plt.xlabel('$k$')
plt.ylabel(r'$|X(k)|$')

plt.subplot(1,2,2)
plt.stem(np.angle(X, True))
plt.grid()
plt.xlabel('$k$')
plt.ylabel(r'$\angle{X(k)}$ (degrees)')

plt.show()
plt.savefig('../figs/X.pdf')
plt.savefig('../figs/X_.eps')
#subprocess.run(shlex.split("termux-open ../figs/X.pdf"))  #if using termex




H = FFT(hn)
plt.figure(2)
plt.stem(hn)
plt.xlabel('$n$')
plt.ylabel('$h(n)$')
plt.grid()
plt.show()
plt.savefig('../figs/h_n.pdf')
plt.savefig('../figs/h_n.eps')
#subprocess.run(shlex.split("termux-open ../figs/h_n.pdf"))  #if using termex

plt.figure(3, figsize=(10,4))

plt.subplot(1,2,1)
plt.stem(np.abs(H))
plt.grid()
plt.xlabel('$k$')
plt.ylabel(r'$|H(k)|$')

plt.subplot(1,2,2)
plt.stem(np.angle(H, True))
plt.grid()
plt.xlabel('$k$')
plt.ylabel(r'$\angle{H(k)}$ (degrees)')

plt.show()
plt.savefig('../figs/H.pdf')
plt.savefig('../figs/H.eps')
#subprocess.run(shlex.split("termux-open ../figs/H.pdf"))  #if using termex


#%%
Y = H * X
yn = IFFT(Y)
plt.figure(4)
plt.stem(np.abs(yn))
plt.xlabel('$n$')
plt.ylabel('$y(n)$')
plt.grid()
plt.show()
plt.savefig('../figs/y_n.pdf')
plt.savefig('../figs/y_n.eps')
#subprocess.run(shlex.split("termux-open ../figs/y_n.pdf"))  #if using termex

plt.figure(5, figsize=(10,4))

plt.subplot(1,2,1)
plt.stem(np.abs(Y))
plt.grid()
plt.xlabel('$k$')
plt.ylabel(r'$|Y(k)|$')

plt.subplot(1,2,2)
plt.stem(np.angle(Y, True))
plt.grid()
plt.xlabel('$k$')
plt.ylabel(r'$\angle{Y(k)}$ (degrees)')

plt.show()
plt.savefig('../figs/Y.pdf')
plt.savefig('../figs/Y.eps')
#subprocess.run(shlex.split("termux-open ../figs/Y.pdf"))  #if using termex


# %%

x = np.random.random(1024)
start_time = timeit.default_timer()
FFT(x)
print(timeit.default_timer() - start_time)

def fft_time():
    SETUP_CODE = ''' 
x = np.random.random(1024)
'''

    TEST_CODE = ''' 
FFT(x)'''

    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    # print mean exec. time
    print('Binary search time: {}'.format(mean(times)))


def DTFT_time():
    SETUP_CODE = ''' 
x = np.random.random(1024)
'''

    TEST_CODE = ''' 
X = DTFT_N2(x)'''

    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    # print mean exec. time
    print('Binary search time: {}'.format(mean(times)))

fft_time()
DTFT_time()
