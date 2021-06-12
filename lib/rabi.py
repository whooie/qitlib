from __future__ import annotations
import numpy as np
import numpy.random as rand
import numpy.fft as fft
import matplotlib.pyplot as pp
import matplotlib.widgets as wig
import lib.plotdefs as pd

def integrator(t: np.ndarray, w0: float, w1: float, w: float,
        noise: float=0) -> tuple[np.ndarray, np.ndarray]:
    """
    Numerically compute the time dependence of the state |psi> = a0*|0> + a1*|1>
    over time grid `t` for relevant frequencies `w0` (atomic frequency), `w1` 
    (Larmor/4*Rabi frequency), and `w2` (drive frequency). Assumes the initial
    population lies entirely in the |0> (ground) state. Optionally adjustable
    `noise` parameter controls the addition of both amplitude and phase noise in
    the drive.
    """
    # add noise with
    #   w1 -> w1 + A \delta_amp
    #   w * t -> w * t + \delta_\phi
    # where \delta_X are normally distributed within some bounds
    dt = abs(t[1] - t[0])
    a0 = np.zeros(t.shape, dtype=np.complex128)
    a1 = np.zeros(t.shape, dtype=np.complex128)
    _w0 = w0 if noise == 0 else w0 * rand.normal(1.0, noise * 1e-3, t.shape)
    _w1 = w1 if noise == 0 else w1 * rand.normal(1.0, noise * 1e-3, t.shape)
    _d = 0 if noise == 0 else rand.normal(0.0, noise * 2*np.pi * 1e-3, t.shape)
    A0 = 2 * dt * (1j * _w1 / 4) * (
        np.exp(1j * (w - _w0) * t + _d) + np.exp(-1j * (w + _w0) * t + _d)
    )
    A1 = 2 * dt * (1j * _w1 / 4) * (
        np.exp(1j * (w + _w0) * t + _d) + np.exp(-1j * (w - _w0) * t + _d)
    )
    a0[:2] = np.cos(w1 / 4 * t[:2])
    a1[:2] = np.sin(w1 / 4 * t[:2]) * 1j
    for i in range(2, t.shape[0]):
        a0[i] = a0[i - 2] + A0[i - 1] * a1[i - 1]
        a1[i] = a1[i - 2] + A1[i - 1] * a0[i - 1]
        if (A := abs(a0[i]) ** 2) + (B := abs(a1[i]) ** 2) > 1:
            a0[i] = a0[i] / np.sqrt(A + B)
            a1[i] = a1[i] / np.sqrt(A + B)
    return a0, a1



