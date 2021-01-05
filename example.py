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

def do_resample():
    fin = 44100
    N = 8192*4
    HN = N // 2
    origSig = np.sin(2.0 * np.pi * 10000.0 / fin * np.arange(N))
    origWin = signal.windows.blackmanharris(N)
    #origWin = np.array([1] * N)
    #origWin /= sum(origWin)
    fftOrig = fft.fft(origSig*origWin, n=N)
    #fftOrig /= sum(origWin)
    fftPlot = np.abs(fftOrig) / sum(origWin)
    freqAxisOrig = fft.fftfreq(N, 1.0 / fin)
#    plt.plot(origSig)
#    plt.show()
#    plt.plot(origWin)
#    plt.show()
#    plt.plot(origSig * origWin)
#    plt.show()
    plt.subplot(211)
    plt.plot(freqAxisOrig[:HN], 20.0 * np.log10(fftPlot)[:HN])
    plt.axis([0, 24000, -200, 0])
    plt.grid(True)
    plt.grid(True, 'both', 'x')
    #plt.show()

    fout = 48000 #192000
    L = 160 #640
    M = 147
    prediv = 8 #2
    coefs = prediv * L * h[::prediv]
    sampleStep = M // L
    output = []
    outputSamples = (N * L) // M
    coefOffset = 0
    sampleIdx = 0
    for i in range(outputSamples):
        cIdx = coefOffset
        sIdx = sampleIdx
        sample = 0.0
        while cIdx < len(coefs):
            sample += coefs[cIdx]*origSig[sIdx]
            cIdx += L
            sIdx -= 1
            if (sIdx < 0): break

        output.append(sample)
        prevOffset = coefOffset
        coefOffset = (coefOffset + M) % L
        sampleIdx = sampleIdx + sampleStep + (1 if prevOffset > coefOffset else 0)

    M = len(output)
    HM = M // 2

    resaSig = np.array(output)
    resaWin = signal.windows.blackmanharris(M)
    fftResa = fft.fft(resaSig*resaWin, n=M)
    fftPlot = np.abs(fftResa) / sum(resaWin)
    freqAxisResa = fft.fftfreq(M, 1.0 / fout)
    plt.subplot(212)
    plt.plot(freqAxisResa[:HM], 20.0 * np.log10(fftPlot)[:HM])
    plt.axis([0, 24000, -200, 0])
    plt.grid(True)
    plt.grid(True, 'both', 'x')
    plt.show()

#plot_filter()
#plot_filter_response()
#plot_coef_table()
do_resample()
