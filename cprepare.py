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

# Compute sinc filter.
h = np.sinc(2 * fL / fS * (np.arange(coefCount) - (coefCount - 1) / 2))
# Apply window.
h = gauss(h)
# Normalize to get unity gain.
#print(np.sum(h))
h /= np.sum(h)

print(h[  1000])
print(h[  8000])
print(h[ 16000])
print(h[ 32000])
print(h[ 64000])
print(h[125000])
print(h[249999])
print(h[250000])
print(h[250001])
print(h[250002])
print(h[250003])
print(h[375002])
print(h[436002])
print(h[468002])
print(h[484002])
print(h[492002])
print(h[499002])

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
