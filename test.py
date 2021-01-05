from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from numpy import fft
from scipy import signal

def log2_up(x):
    result = 1
    while (result < x):
        result <<= 1
    return result

def gauss(x):
    N = len(x)
    HN = (N + 1) // 2
    gaussNum = 0.084 * N
    idx = np.arange(-HN, -HN + N)
    idx = np.divide(idx, gaussNum)
    idx = np.square(idx)
    gaussArr = np.exp(-0.5 * idx)
    return x * gaussArr

# Configuration.
fS = 56448000  # Sampling rate.
fL = 20400  # Cutoff frequency.
coefCount = 500003  # Filter length, must be odd.
beta = 7.857  # Kaiser window beta.

# Compute sinc filter.
h = np.sinc(2 * fL / fS * (np.arange(coefCount) - (coefCount - 1) / 2))

# Apply window.
#h *= np.kaiser(N, beta)
#h *= signal.gaussian(coefCount, std=((1 - 0.6827)*coefCount)/4)
h = gauss(h)

# Normalize to get unity gain.
h /= np.sum(h)

def plot_filter():
    #print(h)
    #print(h[25000:25010])
    plt.plot(h)
    plt.show()

N = log2_up(coefCount)
print(N)

d = h.copy()
d.resize(N)

fftH = fft.fft(d, n=N)
freqaxis = fft.fftfreq(N, 1.0 / fS)

HN = N // 2

def plot_filter_response():
    plt.semilogx()
    plt.plot(freqaxis[:HN], (20.0 * np.log10(np.abs(fftH)))[:HN])
    #plt.plot(freqaxis, (20.0 * np.log10(np.abs(fftH))))
    plt.grid(True)
    plt.grid(True, 'both', 'x')
    plt.show()

def plot_coef_table():
    K = 1280

    coefTable = [[] for _ in range(K)]
    for i, s in enumerate(h):
        coefTable[i % K].append(K * s)

    #print(coefTable[0])

    M = log2_up((N + K - 1) // K)
    freqaxisC = fft.fftfreq(M, 1.0 / (fS / K))
    HM = M // 2
    plt.semilogx()
    for coefs in coefTable:
        fftC = fft.fft(coefs, n=M)
        plt.plot(freqaxisC[:HM], (20.0 * np.log10(np.abs(fftC)))[:HM])
    #    plt.plot(freqaxisC, (20.0 * np.log10(np.abs(fftC))))
    plt.grid(True)
    plt.grid(True, 'both', 'x')
    plt.show()

# Applying the filter to a signal s can be as simple as writing
# s = np.convolve(s, h)

def do_resample(fgen, fin, fout):
    L = fS // fin
    M = fS // fout
    prediv = 1

    while (L % 2) == 0 and (M % 2) == 0:
        prediv *= 2
        L >>= 1
        M >>= 1

    if L < M and (M % L) == 0:
        prediv *= L
        M //= L
        L = 1
    elif M < L and (L % M) == 0:
        prediv *= M
        L //= M
        M = 1
    elif M == L:
        print("No resampling needed")
        return

    print("{}/{}, prediv={}".format(L, M, prediv))

    coefs = prediv * L * h[::prediv]

    axis = [0, max(fin, fout) // 2, -200, 0]

    N = 8192
    HN = N // 2
    outputCount = (N * L) // M
    O = outputCount
    HO = O // 2

    origSig = np.sin(2.0 * np.pi * fgen / fin * np.arange(N))

    sampleStep = M // L
    output = []
    coefOffset = 0
    sampleIdx = 0
    for i in range(outputCount):
        cIdx = coefOffset
        sIdx = sampleIdx
        sample = 0.0 # 0.001 * np.sin(2.0 * np.pi * 1000.0 / fout * i)
        while cIdx < len(coefs):
            sample += coefs[cIdx]*origSig[sIdx]
            cIdx += L
            sIdx -= 1
            if (sIdx < 0):
                break

        output.append(sample)
        prevOffset = coefOffset
        coefOffset = (coefOffset + M) % L
        sampleIdx += sampleStep + (1 if prevOffset > coefOffset else 0)

    resaSig = np.array(output)

    origWin = signal.windows.kaiser(N, beta=33)
    fftOrig = fft.fft(origSig*origWin, n=N)
    origPlot = np.abs(fftOrig / sum(origWin))
    freqAxisOrig = fft.fftfreq(N, 1.0 / fin)

    resaWin = signal.windows.kaiser(O, beta=33)
    fftResa = fft.fft(resaSig*resaWin, n=O)
    resaPlot = np.abs(fftResa / sum(resaWin))
    freqAxisResa = fft.fftfreq(O, 1.0 / fout)

    plt.subplot(211)
    plt.plot(freqAxisOrig[:HN], 20.0 * np.log10(origPlot)[:HN])
    plt.axis(axis)
    plt.grid(True)
    plt.grid(True, 'both', 'x')

    plt.subplot(212)
    plt.plot(freqAxisResa[:HO], 20.0 * np.log10(resaPlot)[:HO])
    plt.axis(axis)
    plt.grid(True)
    plt.grid(True, 'both', 'x')

    plt.show()

#plot_filter()
#plot_filter_response()
#plot_coef_table()
do_resample(11025, 88200, 384000)
