"""
Ref)
    Burrus, Charles & Gopinath, R. & Guo, H.. (1998). Introduction to Wavelets and Wavelet Transform—A Primer. Recherche. 67.
    https://ieeexplore.ieee.org/document/192463
    https://github.com/PyWavelets/pywt/blob/master/pywt/_extensions/c/wavelets_coeffs.template.h
"""

def wavelet(name):
    # default: haar
    leng = 2
    off = 0
    dec_lo = [1 / 2 ** .5, 1 / 2 ** .5]
    dec_hi = [-1 / 2 ** .5, 1 / 2 ** .5]

    if name == 'sym5':
        leng = 10
        off = 0
        dec_lo = [0.019538882735286728, -0.021101834024758855, -0.17532808990845047, 0.016602105764522319,
                  0.63397896345821192, 0.72340769040242059, 0.1993975339773936, -0.039134249302383094,
                  0.029519490925774643, 0.027333068345077982]
        dec_hi = dec_lo[::-1]
        for i in range(0, len(dec_hi), 2):
            dec_hi[i] *= -1

    return leng, off, dec_lo, dec_hi


def dwt_forward(data, wavelet_name='haar', mode='smooth'):
    data = list(data)
    n = len(data)
    if n & 1:
        n ^= 1
        if mode == 'smooth':  # need to compare to pywt
            data.append(data[-1] * 2 - data[-2])
        elif mode == 'zero':
            data.append(0)
        elif mode == 'symmetric':
            data.append(data[-1])
    leng, off, dec_lo, dec_hi = wavelet(wavelet_name)
    low = [0] * (len(data) >> 1)
    high = [0] * (len(data) >> 1)
    nm = leng * n - off
    nh = n >> 1
    for i in range(nh):
        h = g = 0
        k = (i << 1) + nm
        for j in range(leng):
            _k = (k + j) % n
            h += dec_lo[j] * data[_k]
            g += dec_hi[j] * data[_k]
        low[i] += h
        high[i] -= g
    return low, high