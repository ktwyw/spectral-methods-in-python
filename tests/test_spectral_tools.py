"""Unit tests for spectral_tools.  Run with:  pytest tests/"""
import numpy as np
import pytest

import spectral_tools as st


def test_fourier_matrix_matches_fft_and_exact():
    N = 32
    x = st.fourier_grid(N)
    u = np.exp(np.sin(x))
    exact = np.cos(x) * u
    D = st.fourier_diff_matrix(N)
    assert np.abs(D @ u - exact).max() < 1e-12
    assert np.abs(st.fft_derivative(u) - exact).max() < 1e-12
    assert np.abs(D + D.T).max() < 1e-13                      # skew-symmetric


def test_fourier_second_derivative():
    N = 32
    x = st.fourier_grid(N)
    u = np.exp(np.sin(x))
    exact = (np.cos(x) ** 2 - np.sin(x)) * u
    assert np.abs(st.fourier_diff2_matrix(N) @ u - exact).max() < 1e-10
    assert np.abs(st.fft_derivative(u, 2) - exact).max() < 1e-10


def test_cheb_exact_for_polynomials():
    N = 12
    D, x = st.cheb(N)
    for degree in range(N + 1):
        assert np.abs(D @ x**degree - degree * x ** max(degree - 1, 0) * (degree > 0)).max() < 1e-9
    assert np.abs(D.sum(axis=1)).max() < 1e-12


def test_cheb_spectral_accuracy_and_scaling():
    D, x = st.cheb(24)
    u = np.exp(x) * np.sin(3 * x)
    assert np.abs(D @ u - np.exp(x) * (np.sin(3 * x) + 3 * np.cos(3 * x))).max() < 1e-10
    Ds, xs = st.scale_cheb(24, 2.0, 5.0)
    assert xs[0] == pytest.approx(5.0) and xs[-1] == pytest.approx(2.0)
    assert np.abs(Ds @ np.sin(xs) - np.cos(xs)).max() < 1e-10


def test_cheb_coefficients_roundtrip_and_known_series():
    D, x = st.cheb(16)
    a = st.cheb_coeffs(4 * x**3 - 3 * x + 2)                 # = T_3 + 2 T_0
    expected = np.zeros(17)
    expected[0], expected[3] = 2, 1
    assert np.abs(a - expected).max() < 1e-13
    values = np.cos(3 * x) + x
    assert np.abs(st.cheb_values(st.cheb_coeffs(values)) - values).max() < 1e-13


def test_cheb_fft_derivative_matches_matrix():
    D, x = st.cheb(40)
    u = np.exp(-x**2) * np.cos(4 * x)
    assert np.abs(st.cheb_derivative_fft(u) - D @ u).max() < 1e-9


def test_barycentric_interpolation():
    D, x = st.cheb(30)
    f = lambda t: 1 / (1 + 4 * t**2)
    xx = np.linspace(-1, 1, 501)
    assert np.abs(st.barycentric_interpolate(f(x), xx) - f(xx)).max() < 1e-5
    assert np.allclose(st.barycentric_interpolate(f(x), x), f(x))


def test_quadrature_rules():
    x, w = st.clenshaw_curtis(20)
    assert w.sum() == pytest.approx(2.0)
    assert w @ np.exp(x) == pytest.approx(np.e - 1 / np.e, abs=1e-13)
    xg, wg = st.gauss_legendre(10)
    ref_x, ref_w = np.polynomial.legendre.leggauss(10)
    assert np.allclose(xg, ref_x) and np.allclose(wg, ref_w)
    assert wg @ xg**18 == pytest.approx(2 / 19)               # exact to degree 2n - 1 = 19


def test_fft_derivative_nyquist_mode_and_odd_N():
    N = 16
    x = st.fourier_grid(N)
    nyquist = np.cos(N / 2 * x)                               # the sawtooth (-1)^j
    assert np.abs(st.fft_derivative(nyquist, 1)).max() < 1e-12     # odd order: Nyquist mode set to zero
    assert np.abs(st.fft_derivative(nyquist, 2) + (N / 2) ** 2 * nyquist).max() < 1e-10
    assert np.abs(st.fourier_diff2_matrix(N) @ nyquist + (N / 2) ** 2 * nyquist).max() < 1e-10
    x_odd = 2 * np.pi * np.arange(15) / 15                    # odd N has no Nyquist mode
    assert np.abs(st.fft_derivative(np.sin(3 * x_odd)) - 3 * np.cos(3 * x_odd)).max() < 1e-12


def test_fourier_matrices_reject_odd_N():
    for build in (st.fourier_diff_matrix, st.fourier_diff2_matrix):
        with pytest.raises(ValueError):
            build(9)


def test_cheb_small_and_odd_N():
    D0, x0 = st.cheb(0)
    assert D0.shape == (1, 1) and D0[0, 0] == 0 and x0[0] == 1
    D1, x1 = st.cheb(1)
    assert np.allclose(D1, [[0.5, -0.5], [0.5, -0.5]], rtol=0, atol=1e-15)
    D, x = st.cheb(11)                                        # odd N: no grid point at x = 0
    assert np.abs(D @ np.exp(x) - np.exp(x)).max() < 1e-8
    assert np.abs(st.cheb_values(st.cheb_coeffs(np.exp(x))) - np.exp(x)).max() < 1e-13


def test_cheb_coeff_derivative_known_series():
    a = np.zeros(6)
    a[3] = 1.0                                                # T_3' = 3 T_0 + 6 T_2
    assert np.abs(st.cheb_coeff_derivative(a) - np.array([3, 0, 6, 0, 0, 0])).max() < 1e-13


def test_clenshaw_curtis_exact_for_degree_N():
    for N in (1, 2, 7, 8):
        x, w = st.clenshaw_curtis(N)
        for degree in range(N + 1):
            exact = 0.0 if degree % 2 else 2 / (degree + 1)
            assert abs(w @ x**degree - exact) < 1e-13


def test_notebooks_use_the_same_code_as_the_toolbox():
    """Every notebook pastes the functions it needs; the copies must not drift from this module."""
    import ast
    import json
    import pathlib

    def bodies(source):
        found = {}
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.FunctionDef):
                node.returns = None
                for argument in node.args.args:
                    argument.annotation = None
                first = node.body[0]
                if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                    node.body = node.body[1:]                 # ignore docstrings
                found[node.name] = ast.dump(node)
        return found

    root = pathlib.Path(__file__).resolve().parent.parent
    toolbox = bodies((root / "spectral_tools.py").read_text(encoding="utf-8"))
    checked = 0
    for path in sorted(root.glob("notebooks/*.ipynb")) + sorted(root.glob("solutions/*.ipynb")):
        cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
        code = "\n\n".join("".join(cell["source"]) for cell in cells if cell["cell_type"] == "code")
        for name, body in bodies(code).items():
            if name in toolbox:
                assert body == toolbox[name], f"{name} in {path.name} differs from spectral_tools.py"
                checked += 1
    assert checked > 50
