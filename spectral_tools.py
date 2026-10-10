"""Core tools for spectral collocation methods.

A small, dependency-light module (NumPy and SciPy only) that collects the
functions developed step by step in the notebooks of this repository.  Each
notebook defines what it needs inline so that it runs on its own; this file is
the same code in one importable place.

Conventions
-----------
Fourier (periodic) grid :  x_j = 2*pi*j/N,        j = 0, ..., N-1   (N even)
Chebyshev grid          :  x_j = cos(pi*j/N),     j = 0, ..., N     (N+1 points, x_0 = 1, x_N = -1)

Reference for the mathematics: L. N. Trefethen, Spectral Methods in MATLAB,
SIAM, 2000.  The code was written for this repository; `cheb` keeps the name,
the calling convention and the algorithm of the book's cheb.m.

All functions expect one-dimensional real arrays.  The Chebyshev functions that
take values need at least two points (N >= 1) unless stated otherwise.  Note
that fourier_diff2_matrix(N) is not the square of fourier_diff_matrix(N): the
two differ on the highest (Nyquist) mode, as explained in notebook 2.
"""
from __future__ import annotations

import numpy as np
from scipy.fft import dct
from scipy.linalg import eigh_tridiagonal

__all__ = [
    "fourier_grid", "fourier_diff_matrix", "fourier_diff2_matrix", "fft_derivative",
    "cheb", "cheb_coeffs", "cheb_values", "cheb_coeff_derivative", "cheb_derivative_fft",
    "barycentric_interpolate", "clenshaw_curtis", "gauss_legendre", "scale_cheb",
]


# ----------------------------------------------------------------------------
# Periodic problems: Fourier
# ----------------------------------------------------------------------------
def fourier_grid(N: int) -> np.ndarray:
    """Equispaced periodic grid on [0, 2*pi): x_j = 2*pi*j/N, j = 0..N-1."""
    return 2 * np.pi * np.arange(N) / N


def fourier_diff_matrix(N: int) -> np.ndarray:
    """First-derivative matrix for the periodic grid of N points (N even).

    Entry (i, j) is (1/2) (-1)^(i-j) cot((x_i - x_j)/2) off the diagonal and 0 on it.
    """
    if N % 2:
        raise ValueError("N must be even")
    h = 2 * np.pi / N
    k = np.arange(N)[:, None] - np.arange(N)[None, :]          # i - j
    D = np.zeros((N, N))
    off = k != 0
    D[off] = 0.5 * (-1.0) ** k[off] / np.tan(k[off] * h / 2)
    return D


def fourier_diff2_matrix(N: int) -> np.ndarray:
    """Second-derivative matrix for the periodic grid of N points (N even)."""
    if N % 2:
        raise ValueError("N must be even")
    h = 2 * np.pi / N
    k = np.arange(N)[:, None] - np.arange(N)[None, :]
    D2 = np.full((N, N), -np.pi**2 / (3 * h**2) - 1 / 6)
    off = k != 0
    D2[off] = -0.5 * (-1.0) ** k[off] / np.sin(k[off] * h / 2) ** 2
    return D2


def fft_derivative(v: np.ndarray, order: int = 1) -> np.ndarray:
    """Spectral derivative of real periodic samples v on fourier_grid(len(v)).

    For odd orders the highest (Nyquist) mode is set to zero, because on the
    grid it cannot be told apart from its mirror image and its derivative is
    ambiguous.
    """
    N = len(v)
    k = np.fft.fftfreq(N, d=1 / N)                               # 0, 1, ..., N/2-1, -N/2, ..., -1
    factor = (1j * k) ** order
    if order % 2 and N % 2 == 0:
        factor[N // 2] = 0
    return np.real(np.fft.ifft(factor * np.fft.fft(v)))


# ----------------------------------------------------------------------------
# Bounded intervals: Chebyshev
# ----------------------------------------------------------------------------
def cheb(N: int) -> tuple[np.ndarray, np.ndarray]:
    """Chebyshev differentiation matrix D ((N+1) x (N+1)) and points x_j = cos(pi*j/N).

    Off-diagonal entries are (c_i/c_j) (-1)^(i+j) / (x_i - x_j) with c_0 = c_N = 2
    and c_j = 1 otherwise.  The diagonal is set so that every row sums to zero
    (the derivative of a constant is zero), which is more accurate in floating
    point than the explicit diagonal formulas.
    """
    if N == 0:
        return np.zeros((1, 1)), np.ones(1)
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    c = np.where((j == 0) | (j == N), 2.0, 1.0) * (-1.0) ** j
    dx = x[:, None] - x[None, :]
    D = (c[:, None] / c[None, :]) / (dx + np.eye(N + 1))         # the identity avoids 0/0 on the diagonal
    D -= np.diag(D.sum(axis=1))
    return D, x


def scale_cheb(N: int, a: float, b: float) -> tuple[np.ndarray, np.ndarray]:
    """Chebyshev matrix and points mapped from [-1, 1] to [a, b] (x runs from b down to a)."""
    D, x = cheb(N)
    return D * 2 / (b - a), a + (b - a) * (x + 1) / 2


def cheb_coeffs(values: np.ndarray) -> np.ndarray:
    """Chebyshev coefficients a_k of the degree-N interpolant through values at cheb points.

    p(x) = sum_k a_k T_k(x).  Computed with a type-I discrete cosine transform.
    """
    values = np.asarray(values, dtype=float)
    N = len(values) - 1
    if N == 0:
        return values.copy()
    a = dct(values, type=1) / N
    a[0] /= 2
    a[N] /= 2
    return a


def cheb_values(coeffs: np.ndarray) -> np.ndarray:
    """Inverse of cheb_coeffs: values at the Chebyshev points from coefficients."""
    a = np.asarray(coeffs, dtype=float).copy()
    N = len(a) - 1
    if N == 0:
        return a
    a[0] *= 2
    a[N] *= 2
    return dct(a, type=1) / 2


def cheb_coeff_derivative(a: np.ndarray) -> np.ndarray:
    """Coefficients of p'(x) from the coefficients of p(x), by the backward recurrence."""
    a = np.asarray(a, dtype=float)
    N = len(a) - 1
    b = np.zeros(N + 1)
    if N == 0:
        return b
    b[N - 1] = 2 * N * a[N]
    for k in range(N - 2, -1, -1):
        b[k] = b[k + 2] + 2 * (k + 1) * a[k + 1]
    b[0] /= 2
    return b


def cheb_derivative_fft(values: np.ndarray) -> np.ndarray:
    """Derivative at the Chebyshev points in O(N log N): transform, differentiate, transform back."""
    return cheb_values(cheb_coeff_derivative(cheb_coeffs(values)))


def barycentric_interpolate(values: np.ndarray, x_eval: np.ndarray) -> np.ndarray:
    """Evaluate the polynomial interpolant through values at the Chebyshev points.

    Uses the barycentric formula, which is numerically stable for any N.
    """
    values = np.asarray(values, dtype=float)
    N = len(values) - 1
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    w = (-1.0) ** j
    w[0] /= 2
    w[N] /= 2
    x_eval = np.atleast_1d(np.asarray(x_eval, dtype=float))
    diff = x_eval[:, None] - x[None, :]
    exact = diff == 0
    diff[exact] = 1.0
    ratio = w / diff
    result = ratio @ values / ratio.sum(axis=1)
    hit_rows, hit_cols = np.nonzero(exact)
    result[hit_rows] = values[hit_cols]
    return result


# ----------------------------------------------------------------------------
# Integration
# ----------------------------------------------------------------------------
def clenshaw_curtis(N: int) -> tuple[np.ndarray, np.ndarray]:
    """Clenshaw-Curtis nodes (the Chebyshev points) and weights on [-1, 1].

    The weights integrate the degree-N interpolant exactly.  They follow from
    the integrals of the Chebyshev polynomials: int T_k = 2/(1 - k^2) for even k.
    """
    if N == 0:
        return np.ones(1), np.array([2.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    moments = np.zeros(N + 1)
    even = np.arange(0, N + 1, 2)
    moments[even] = 2 / (1 - even.astype(float) ** 2)
    # integral of the interpolant = sum_k a_k * moments_k, and a = cheb_coeffs(values) is a
    # cosine transform of the values; transposing that transform gives the weights
    w = dct(moments, type=1) / N
    w[0] /= 2
    w[N] /= 2
    return x, w


def gauss_legendre(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes and weights (n points) by the Golub-Welsch eigenvalue method."""
    k = np.arange(1, n)
    beta = k / np.sqrt(4.0 * k**2 - 1)
    nodes, vectors = eigh_tridiagonal(np.zeros(n), beta)
    return nodes, 2 * vectors[0] ** 2
