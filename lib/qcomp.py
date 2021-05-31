from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as pp
import lib.plotdefs as pd
from lib.dirac import *
import re
import copy
import random
import pathlib

"""
Provides machinery to construct and run quantum computing circuits using objects
from lib.dirac. Supported gates are:
    NOT, CNOT, PHASE, CPHASE, HADAMARD, SWAP, REVERSE, QFT, IQFT
and are implemented all in terms of the H/P/CNOT universal gate set. Gates
    XYMODN, CXYMODN
are also supported for Shor's algorithm, but not implemented as H/P/CNOT gates.

Also provides simple functions for converting between binary and decimal
representations of numbers. Bitlists are in big-endian order, while bitstrs and
qubit state labels are in small-endian order.
"""
### bitlists are in big-endian order, while bitstrs are in small-endian order

def int_to_bitstr(k: int, n: int=None):
    b = str(bin(k))[2:]
    _n = len(b) if n is None or n < 0 else n
    return ((_n - len(b)) * "0" + b)[-_n:]

def i2b(k: int, n: int):
    return int_to_bitstr(k, n)

def bitstr_to_int(B: str):
    return int(B, 2)

def b2i(B: str):
    return bitstr_to_int(B)

def qubit_to_int(S: BasisState):
    return bitstr_to_int(S.label)

def int_to_bitlist(k: int, n: int=None):
    return [int(b) for b in int_to_bitstr(k, n)][::-1]

def bitlist_to_bitstr(B: list):
    return "".join(str(b) for b in B[::-1])

def bitlist_to_int(B: list):
    return bitstr_to_int(bitlist_to_bitstr(B))

def qubit_basis(n: int):
    """
    Generate a `Basis` for `n` qubit wires.
    """
    return Basis(*[
        BasisState(decstr_to_bitstr(str(k), n)) for k in range(2**n)
    ])

def _hadamard_vecop(n: int, k: int):
    _action = dict()
    for i in range(2**n):
        I = i2b(i, n)
        s0 = BasisState(I[:k] + "0" + I[k + 1:])
        s1 = BasisState(I[:k] + "1" + I[k + 1:])
        _action.update({
            s0: StateVec({s0: +1/np.sqrt(2), s1: +1/np.sqrt(2)}, is_ket=True),
            s1: StateVec({s0: +1/np.sqrt(2), s1: -1/np.sqrt(2)}, is_ket=True)
        })
    return VecOperator(_action, default_zero=False, is_ketop=True, scalar=1.0)

def _hadamard_mat(n: int, k: int, sparse: bool=False):
    if n is None:
        raise ValueError("Wire count cannot be `None`")
    H = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    acc = sp._csr_matrix([[1]])
    for j in range(n):
        acc = sp.kron(acc, H if j == k else sp.eye(2), format="csr")
    return acc if sparse else acc.toarray()

def _hadamard_matop(n: int, k: int, sparse: bool=False):
    return MatOperator(
        _hadamard_mat(n, k, sparse),
        qubit_basis(n),
        is_ketop=True
    )

def _hadamard_func(b: BasisState, k: int):
    l = b.label
    if l[k] == "0":
        return StateVec({
            BasisState(l[:k]+"0"+l[k+1:]): +1/np.sqrt(2),
            BasisState(l[:k]+"1"+l[k+1:]): +1/np.sqrt(2)
        }, is_ket=True)
    elif b.label[k] == "1":
        return StateVec({
            BasisState(l[:k]+"0"+l[k+1:]): +1/np.sqrt(2),
            BasisState(l[:k]+"1"+l[k+1:]): -1/np.sqrt(2)
        }, is_ket=True)

def _hadamard_funcop(n: int, k: int):
    return FuncOperator(
        lambda b: _hadamard_func(b, k),
        is_ketop=True
    )

def hadamard(n: int, k: int, kind: str="mat", **kwargs):
    """
    Construct a `kind`-`Operator` representing a Hadamard gate operating on the
    `k`-th qubit (out of `n` total qubits).
    """
    if kind == "vec":
        return _hadamard_vecop(n, k)
    elif kind == "mat":
        return _hadamard_matop(n, k, kwargs.get("sparse", False))
    elif kind == "func":
        return _hadamard_funcop(n, k)

def _phase_vecop(n: int, k: int, theta: float):
    _action = dict()
    for i in range(2**n):
        I = i2b(i, n)
        s0 = BasisState(I[:k] + "0" + I[k + 1:])
        s1 = BasisState(I[:k] + "1" + I[k + 1:])
        _action.update({
            s0: StateVec({s0: 1.0 + 0.0j}, is_ket=True),
            s1: StateVec({s1: np.exp(1.0j * theta)}, is_ket=True)
        })
    return VecOperator(_action, default_zero=False, is_ketop=True, scalar=1.0)

def _phase_mat(n: int, k: int, theta: float, sparse: bool=False):
    if n is None:
        raise ValueError("Wire count cannot be `None`")
    P = np.array([[1, 0], [0, np.exp(1j*theta)]])
    acc = sp._csr_matrix([[1]])
    for j in range(n):
        acc = sp.kron(acc, P if j == k else sp.eye(2), format="csr")
    return acc if sparse else acc.toarray()

def _phase_matop(n: int, k: int, theta: float, sparse: bool=False):
    return MatOperator(
        _phase_mat(n, k, theta, sparse),
        qubit_basis(n),
        is_ketop=True
    )

def _phase_func(b: BasisState, k: int, theta: float):
    l = b.label
    if l[k] == "1":
        return StateVec({b: np.exp(1.0j * theta)}, is_ket=True)
    else:
        return StateVec({b: 1.0 + 0.0j}, is_ket=True)

def _phase_funcop(n: int, k: int, theta: float):
    return FuncOperator(
        lambda b: _phase_func(b, k, theta),
        is_ketop=True
    )

def phase(n: int, k: int, theta: float, kind: str="mat", **kwargs):
    """
    Construct a `kind`-`Operator` representing a `theta`-phase gate operating on
    the `k`-th qubit (out of `n` total qubits).
    """
    if kind == "vec":
        return _phase_vecop(n, k, theta)
    elif kind == "mat":
        return _phase_matop(n, k, theta, kwargs.get("sparse", False))
    elif kind == "func":
        return _phase_funcop(n, k, theta)

def _cnot_vecop(n: int, k1: int, k2: int):
    assert k1 != k2
    j1 = min(k1, k2)
    j2 = max(k1, k2)
    b = ("1", "0") if k2 > k1 else ("0", "1")
    B = ("1", "1")
    _action = dict()
    for i in range(2**n):
        I = i2b(i, n)
        s1 = BasisState(I[:j1] + b[0] + I[j1 + 1:j2] + b[1] + I[j2 + 1:])
        s2 = BasisState(I[:j1] + B[0] + I[j1 + 1:j2] + B[1] + I[j2 + 1:])
        _action.update({s1: s2, s2: s1})
    return VecOperator(_action, default_zero=False, is_ketop=True, scalar=1.0)

def _cnot_mat(n: int, k1: int, k2: int, sparse: bool=False):
    assert k1 != k2
    if n is None:
        raise ValueError("Wire count cannot be `None`")
    X = sp._csr_matrix((2**n, 2**n), dtype=complex) if sparse \
            else np.zeros((2**n, 2**n), dtype=complex)
    for k in range(2**n):
        B = int_to_bitlist(k, n)
        B[k2] = (B[k1] + B[k2]) % 2
        b = bitlist_to_int(B)
        X[b, k] = 1
    return X

def _cnot_matop(n: int, k1: int, k2: int, sparse: bool=False):
    return MatOperator(
        _cnot_mat(n, k1, k2, sparse),
        qubit_basis(n),
        is_ketop=True
    )

def _cnot_func(b: BasisState, k1: int, k2: int):
    l = b.label
    if l[k1] == "1":
        return StateVec({
            BasisState(l[:k2]+("0" if l[k2] == "1" else "1")+l[k2+1:]): \
                    1.0 + 0.0j
        })
    else:
        return StateVec({b: 1.0 + 0.0j})

def _cnot_funcop(n: int, k1: int, k2: int):
    return FuncOperator(
        lambda b: _cnot_func(b, k1, k2),
        is_ketop=True
    )

def cnot(n: int, k1: int, k2: int, kind: str="mat", **kwargs):
    """
    Construct a `kind`-`Operator` representing a CNOT gate operating on the
    `k1`-th qubit with the `k2`-th qubit as control (out of `n` total qubits).
    """
    if kind == "vec":
        return _cnot_vecop(n, k1, k2)
    elif kind == "mat":
        return _cnot_matop(n, k1, k2, kwargs.get("sparse", False))
    elif kind == "func":
        return _cnot_funcop(n, k1, k2)

def _qft_func(b: BasisState):
    n = len(b.label)
    N = 2**n
    j = bitstr_to_int(b.label)
    return StateVec({
        BasisState(i2b(k, n)): np.exp(2*np.pi*1j * j*k/N) / np.sqrt(N) \
                for k in range(N)
    })

def qft_elements(k1: int, k2: int):
    """
    Generate the circuit elements for the QFT acting on all wires from `k1` to
    `k2` (inclusive). Assumes `k2` >= `k1`.
    """
    assert k1 <= k2
    if k1 == k2:
        return [("H", k2)]
    else:
        elements = qft_elements(k1 + 1, k2) \
                + [("CP", k2 - i, k1, np.pi / 2**(k2 - k1 - i))
                    for i in range(k2 - k1)] \
                + [("H", k1)]
        return elements

def iqft_elements(k1: int, k2: int):
    """
    Generate the circuit elements for the inverse QFT acting on all wires from
    `k1` to `k2` (inclusive). Assumes `k2` >= `k1`.
    """
    assert k1 <= k2
    if k1 == k2:
        return [("H", k2)]
    else:
        elements = [("H", k1)] \
                + [("CP", k2 - i, k1, -np.pi / 2**(k2 - k1 - i))
                    for i in reversed(range(k2 - k1))] \
                + iqft_elements(k1 + 1, k2)
        return elements

def reverse_elements(k1: int, k2: int):
    """
    Generate the circuit elements to reverse the order of all wires from `k1` to
    `k2` (inclusive). Assumes `k2` >= `k1`.
    """
    assert k1 <= k2
    elements = list()
    for i in range((k2 - k1 + 1) // 2):
        elements += [
            ("CNOT", k1 + i, k2 - i),
            ("CNOT", k2 - i, k1 + i),
            ("CNOT", k1 + i, k2 - i)
        ]
    return elements

def _xymodN_vecop(n: int, k: int, q: int, x: int, N: int):
    assert 2**q >= N
    assert k + q - 1 < n
    _action = dict()
    for i in range(2**n):
        I = i2b(i, n)
        j = b2i(I[k:k + q])
        s1 = BasisState(I)
        if j < N:
            s2 = BasisState(I[:k] + i2b((j * x) % N, q) + I[k + q:])
        else:
            s2 = s1
        _action[s1] = StateVec({s2: 1.0 + 0.0j}, is_ket=True)
    return VecOperator(_action, default_zero=False, is_ketop=True, scalar=1.0)

def _xymodN_mat(n: int, k: int, q: int, x: int, N: int, sparse: bool=False):
    assert 2**q >= N
    assert k + q - 1 < n
    X = sp._csr_matrix((2**n, 2**n), dtype=complex) if sparse \
            else np.zeros((2**n, 2**n), dtype=complex)
    for i in range(2**n):
        I = i2b(i, n)
        j = b2i(I[k:k + q])
        if j < N:
            J = (j * x) % N
        else:
            J = j
        X[b2i(I[:k] + i2b(J, q) + I[k + q:]), i] = 1.0 + 0.0j
    return X

def _xymodN_matop(n: int, k: int, q: int, x: int, N: int, sparse: bool=False):
    return MatOperator(
        _xymodN_mat(n, k, q, x, N, sparse),
        qubit_basis(n),
        is_ketop=True
    )

def _xymodN_func(b: BasisState, k: int, q: int, x: int, N: int):
    assert 2**q >= N
    I = b.label
    j = b2i(I[k:k + q])
    if j < N:
        J = (j * x) % N
        return StateVec(
            {BasisState(I[:k] + i2b(J, q) + I[k + q:]): 1.0 + 0.0j},
            is_ket=True
        )
    else:
        return StateVec({b: 1.0 + 0.0j}, is_ket=True)

def _xymodN_funcop(n: int, k: int, q: int, x: int, N: int):
    assert k + q - 1 < n
    return FuncOperator(
        lambda b: _xymodN_func(b, k, q, x, N),
        is_ketop=True
    )

def xymodN(n: int, k: int, q: int, x: int, N: int, kind: str="mat", **kwargs):
    """
    Construct a `kind`-`Operator` representing a gate operating on qubits `k`
    through `k + q - 1` (out of `n` total qubits), performing the operation
        x * y mod N
    where `x` and `N` are fixed constants and `y` is the numerical value of the
    specified qubits.
    """
    if kind == "vec":
        return _xymodN_vecop(n, k, q, x, N)
    elif kind == "mat":
        return _xymodN_matop(n, k, q, x, N, kwargs.get("sparse", False))
    elif kind == "func":
        return _xymodN_funcop(n, k, q, x, N)

def _cxymodN_vecop(n: int, k1: int, k2: int, q: int, x: int, N: int):
    assert 2**q >= N
    assert k2 + q - 1 < n
    assert k1 < k2 or k1 >= k2 + q
    _action = dict()
    for i in range(2**n):
        I = i2b(i, n)
        j = b2i(I[k2:k2 + q])
        s1 = BasisState(I)
        if j < N and I[k1] == "1":
            s2 = BasisState(I[:k2] + i2b((j * x) % N, q) + I[k2 + q:])
        else:
            s2 = s1
        _action[s1] = StateVec({s2: 1.0 + 0.0j}, is_ket=True)
    return VecOperator(_action, default_zero=False, is_ketop=True, scalar=1.0)

def _cxymodN_mat(n: int, k1: int, k2: int, q: int, x: int, N: int,
        sparse: bool=False):
    assert 2**q >= N
    assert k2 + q - 1 < n
    assert k1 < k2 or k1 >= k2 + q
    X = sp._csr_matrix((2**n, 2**n), dtype=complex) if sparse \
            else np.zeros((2**n, 2**n), dtype=complex)
    for i in range(2**n):
        I = i2b(i, n)
        j = b2i(I[k2:k2 + q])
        if j < N and I[k1] == "1":
            J = (j * x) % N
        else:
            J = j
        X[b2i(I[:k2] + i2b(J, q) + I[k2 + q:]), i] = 1.0 + 0.0j
    return X

def _cxymodN_matop(n: int, k1: int, k2: int, q: int, x: int, N: int,
        sparse: bool=False):
    return MatOperator(
        _cxymodN_mat(n, k1, k2, q, x, N, sparse),
        qubit_basis(n),
        is_ketop=True
    )

def _cxymodN_func(b: BasisState, k1: int, k2: int, q: int, x: int, N: int):
    assert 2**q >= N
    assert k1 < k2 or k1 >= k2 + q
    I = b.label
    j = b2i(I[k2:k2 + q])
    if j < N and I[k1] == "1":
        J = (j * x) % N
        return StateVec(
            {BasisState(I[:k2] + i2b(J, q) + I[k2 + q:]): 1.0 + 0.0j},
            is_ket=True
        )
    else:
        return StateVec({b: 1.0 + 0.0j}, is_ket=True)

def _cxymodN_funcop(n: int, k1: int, k2: int, q: int, x: int, N: int):
    assert k2 + q - 1 < n
    return FuncOperator(
        lambda b: _cxymodN_func(b, k1, k2, q, x, N),
        is_ketop=True
    )

def cxymodN(n: int, k1: int, k2: int, q: int, x: int, N: int, kind: str="mat",
        **kwargs):
    """
    Construct a `kind`-`Operator` representing a gate operating on qubits `k1`
    through `k1 + q - 1` (out of `n` total qubits) with qubit `k2` as the
    control, performing the operation
        x * y mod N
    where `x` and `N` are fixed constants and `y` is the numerical value of the
    specified qubits.
    """
    if kind == "vec":
        return _cxymodN_vecop(n, k1, k2, q, x, N)
    elif kind == "mat":
        return _cxymodN_matop(n, k1, k2, q, x, N, kwargs.get("sparse", False))
    elif kind == "func":
        return _cxymodN_funcop(n, k1, k2, q, x, N)

class CircuitOutput:
    """
    Stores and manages the output state and measurement of a quantum circuit.
    """
    def __init__(self, state: StateVec | StateArr, basis: Basis,
            measurement: list[int]=None):
        self.state = state
        self.basis = basis
        self.measurement = measurement

    def get_statevec(self):
        return self.state if isinstance(self.state, StateVec) \
                else self.state.to_statevec()

    def get_statearr(self):
        return self.state if isinstance(self.state, StateArr) \
                else self.state.to_statearr(self.basis)

    def get_measurement(self):
        return self.measurement

    def get_hist(self, density: bool=False, fractional: bool=False):
        if self.measurement is None:
            return None
        bins = np.arange(-0.5, len(self.basis) + 0.5, 1)
        if fractional:
            measurement = list(map(
                lambda x: x / len(self.basis),
                self.measurement
            ))
            bins = bins / len(self.basis)
        else:
            measurement = self.measurement
        P = pd.Plotter()
        P \
            .hist(measurement, bins=bins, edgecolor="k", density=density) \
            .set_xlabel("Measurement") \
            .set_ylabel("Probability density" if density else "Count") \
            .grid(False, "both", "x")
        return P

    def save_hist(self, outfilename: str, density: bool=False,
            fractional: bool=False):
        if self.measurement is None:
            raise Exception("Does not contain measurement")
        else:
            P = self.get_hist(density, fractional)
            P.savefig(outfilename)

class Circuit:
    """
    Drives the creation and simulation of a quantum circuit. The circuit itself
    is stored as a simple list of gate names and associated arguments, rather
    than as full `*Operator` objects to save on memory. This list can be
    compiled to a list of only gates that are implemented atomically or to a
    single `*Operator` representing the overall action of the entire circuit.

    `Circuit`s may be created from a text file detailing the initial state and
    various gates of the circuit, and may also be used to write similar text
    files.

    When run, `Circuit`s may perform gate operations one-by-one or by first
    generating the overall operator mentioned above using any of the kinds of
    `*Operator`s defined in lib.dirac. `MatOperator`s may be used with sparse
    internal representations as well.

    If a measurement at the end of a `Circuit` is desired, a number of
    measurements to perform may optionally be specified. If a measurement is
    desired but this number is left unspecified, the default number is
        2^N * 1000
    where `N` is the number of qubits used by the `Circuit`.
    """
    def __init__(self, nwires: int, elements: list[tuple],
            initial_state: StateArr, measure: bool | int=False):
        self.nwires = nwires
        self.basis = qubit_basis(nwires)
        self.elements = elements
        self.initial_state = initial_state
        self.measure = measure

    def __getitem__(self, pos: int):
        return self.elements[pos]

    @staticmethod
    def get_vecoperator(gate: str, nwires: int, args: tuple):
        """
        Generate the appropriate `VecOperator` for a given gate name and
        associated arguments.
        """
        if gate in {"H", "HADAMARD"} and len(args) == 1:
            return _hadamard_vecop(nwires, *args)
        elif gate in {"P", "PHASE"} and len(args) == 2:
            return _phase_vecop(nwires, *args)
        elif gate == "CNOT" and len(args) == 2:
            return _cnot_vecop(nwires, *args)
        elif gate == "NOT" and len(args) == 1:
            (k,) = args
            return _hadamard_vecop(nwires, k) \
                    * _phase_vecop(nwires, k, np.pi) \
                    * _hadamard_vecop(nwires, k)
        elif gate == "RZ" and len(args) == 2:
            (k, theta) = args
            return _hadamard_vecop(nwires, k) \
                    * _phase_vecop(nwires, k, np.pi) \
                    * _hadamard_vecop(nwires, k) \
                    * _phase_vecop(nwires, k, -theta/2) \
                    * _hadamard_vecop(nwires, k) \
                    * _phase_vecop(nwires, k, np.pi) \
                    * _hadamard_vecop(nwires, k) \
                    * _phase_vecop(nwires, k, theta/2)
        elif gate == "CRZ" and len(args) == 3:
            (k1, k2, theta) = args
            return _cnot_vecop(nwires, k1, k2) \
                    * _phase_vecop(nwires, k2, -theta/2) \
                    * _cnot_vecop(nwires, k1, k2) \
                    * _phase_vecop(nwires, k2, theta/2)
        elif gate in {"CP", "CPHASE"} and len(args) == 3:
            (k1, k2, theta) = args
            return _phase_vecop(nwires, k1, theta/2) \
                    * _cnot_vecop(nwires, k1, k2) \
                    * _phase_vecop(nwires, k2, -theta/2) \
                    * _cnot_vecop(nwires, k1, k2) \
                    * _phase_vecop(nwires, k2, theta/2)
        elif gate == "SWAP" and len(args) == 2:
            (k1, k2) = args
            return _cnot_vecop(nwires, k1, k2) \
                    * _cnot_vecop(nwires, k2, k1) \
                    * _cnot_vecop(nwires, k1, k2)
        elif gate == "QFT" and len(args) == 2:
            (k1, k2) = args
            acc = VecOperator.identity(is_ketop=True)
            for X in qft_elements(k1, k2):
                acc = Circuit.get_vecoperator(X[0], nwires, X[1:]) * acc
            return acc
        elif gate == "IQFT" and len(args) == 2:
            (k1, k2) = args
            acc = VecOperator.identity(is_ketop=True)
            for X in iqft_elements(k1, k2):
                acc = Circuit.get_vecoperator(X[0], nwires, X[1:]) * acc
            return acc
        elif gate == "REVERSE" and len(args) == 2:
            (k1, k2) = args
            acc = VecOperator.identity(is_ketop=True)
            for X in reverse_elements(k1, k2):
                acc = Circuit.get_vecoperator(X[0], nwires, X[1:]) * acc
            return acc
        elif gate == "FUNC" and len(args) == 5 and args[2] == "xyModN":
            (k, q, _, x, N) = args
            return _xymodN_vecop(nwires, k, q, x, N)
        elif gate == "CFUNC" and len(args) == 6 and args[3] == "xyModN":
            (k1, k2, q, _, x, N) = args
            return _cxymodN_vecop(nwires, k1, k2, q, x, N)
        else:
            raise Exception("Invalid circuit element")

    @staticmethod
    def get_mat(gate: str, nwires: int, args: tuple, sparse: bool=False):
        """
        Generate the appropriate matrix for a given gate name and associated
        arguments.
        """
        if gate in {"H", "HADAMARD"} and len(args) == 1:
            return _hadamard_mat(nwires, *args, sparse)
        elif gate in {"P", "PHASE"} and len(args) == 2:
            return _phase_mat(nwires, *args, sparse)
        elif gate == "CNOT" and len(args) == 2:
            return _cnot_mat(nwires, *args, sparse)
        elif gate == "NOT" and len(args) == 1:
            (k,) = args
            return _hadamard_mat(nwires, k, sparse) \
                    @ _phase_mat(nwires, k, np.pi, sparse) \
                    @ _hadamard_mat(nwires, k, sparse)
        elif gate == "RZ" and len(args) == 2:
            (k, theta) = args
            return _hadamard_mat(nwires, k, sparse) \
                    @ _phase_mat(nwires, k, np.pi, sparse) \
                    @ _hadamard_mat(nwires, k, sparse) \
                    @ _phase_mat(nwires, k, -theta/2, sparse) \
                    @ _hadamard_mat(nwires, k, sparse) \
                    @ _phase_mat(nwires, k, np.pi, sparse) \
                    @ _hadamard_mat(nwires, k, sparse) \
                    @ _phase_mat(nwires, k, theta/2, sparse)
        elif gate == "CRZ" and len(args) == 3:
            (k1, k2, theta) = args
            return _cnot_mat(nwires, k1, k2, sparse) \
                    @ _phase_mat(nwires, k2, -theta/2, sparse) \
                    @ _cnot_mat(nwires, k1, k2, sparse) \
                    @ _phase_mat(nwires, k2, theta/2, sparse)
        elif gate in {"CP", "CPHASE"} and len(args) == 3:
            (k1, k2, theta) = args
            return _phase_mat(nwires, k1, theta/2, sparse) \
                    @ _cnot_mat(nwires, k1, k2, sparse) \
                    @ _phase_mat(nwires, k2, -theta/2, sparse) \
                    @ _cnot_mat(nwires, k1, k2, sparse) \
                    @ _phase_mat(nwires, k2, theta/2, sparse)
        elif gate == "SWAP" and len(args) == 2:
            (k1, k2) = args
            return _cnot_mat(nwires, k1, k2, sparse) \
                    @ _cnot_mat(nwires, k2, k1, sparse) \
                    @ _cnot_mat(nwires, k1, k2, sparse)
        elif gate == "QFT" and len(args) == 2:
            (k1, k2) = args
            acc = (sp.identity(2**nwires, dtype=complex, format="csr") \
                    if sparse else np.identity(2**nwires, dtype=complex))
            for X in qft_elements(k1, k2):
                acc = Circuit.get_mat(X[0], nwires, X[1:], sparse) @ acc
            return acc
        elif gate == "IQFT" and len(args) == 2:
            (k1, k2) = args
            acc = (sp.identity(2**nwires, dtype=complex, format="csr") \
                    if sparse else np.identity(2**nwires, dtype=complex))
            for X in iqft_elements(k1, k2):
                acc = Circuit.get_mat(X[0], nwires, X[1:], sparse) @ acc
            return acc
        elif gate == "REVERSE" and len(args) == 2:
            (k1, k2) = args
            acc = (sp.identity(2**nwires, dtype=complex, format="csr") \
                    if sparse else np.identity(2**nwires, dtype=complex))
            for X in reverse_elements(k1, k2):
                acc = Circuit.get_mat(X[0], nwires, X[1:], sparse) @ acc
            return acc
        elif gate == "FUNC" and len(args) == 5 and args[2] == "xyModN":
            (k, q, _, x, N) = args
            return _xymodN_mat(nwires, k, q, x, N)
        elif gate == "CFUNC" and len(args) == 6 and args[3] == "xyModN":
            (k1, k2, q, _, x, N) = args
            return _cxymodN_mat(nwires, k1, k2, q, x, N)
        else:
            raise Exception("Invalid circuit element")

    @staticmethod
    def get_matoperator(gate: str, nwires: int, args: tuple,
            sparse: bool=False):
        """
        Generate the appropriate `MatOperator` for a given gate name and
        associated arguments.
        """
        return MatOperator(
            Circuit.get_mat(gate, nwires, args, sparse),
            qubit_basis(nwires),
            is_ketop=True
        )

    @staticmethod
    def get_funcoperator(gate: str, nwires: int, args: tuple):
        """
        Generate the appropriate `FuncOperator` for a given gate name and
        associated arguments.
        """
        if gate in {"H", "HADAMARD"} and len(args) == 1:
            return _hadamard_funcop(nwires, *args)
        elif gate in {"P", "PHASE"} and len(args) == 2:
            return _phase_funcop(nwires, *args)
        elif gate == "CNOT" and len(args) == 2:
            return _cnot_funcop(nwires, *args)
        elif gate == "NOT" and len(args) == 1:
            (k,) = args
            return _hadamard_funcop(nwires, k) \
                    * _phase_funcop(nwires, k, np.pi) \
                    * _hadamard_funcop(nwires, k)
        elif gate == "RZ" and len(args) == 2:
            (k, theta) = args
            return _hadamard_funcop(nwires, k) \
                    * _phase_funcop(nwires, k, np.pi) \
                    * _hadamard_funcop(nwires, k) \
                    * _phase_funcop(nwires, k, -theta/2) \
                    * _hadamard_funcop(nwires, k) \
                    * _phase_funcop(nwires, k, np.pi) \
                    * _hadamard_funcop(nwires, k) \
                    * _phase_funcop(nwires, k, theta/2)
        elif gate == "CRZ" and len(args) == 3:
            (k1, k2, theta) = args
            return _cnot_funcop(nwires, k1, k2) \
                    * _phase_funcop(nwires, k2, -theta/2) \
                    * _cnot_funcop(nwires, k1, k2) \
                    * _phase_funcop(nwires, k2, theta/2)
        elif gate in {"CP", "CPHASE"} and len(args) == 3:
            (k1, k2, theta) = args
            return _phase_funcop(nwires, k1, theta/2) \
                    * _cnot_funcop(nwires, k1, k2) \
                    * _phase_funcop(nwires, k2, -theta/2) \
                    * _cnot_funcop(nwires, k1, k2) \
                    * _phase_funcop(nwires, k2, theta/2)
        elif gate == "SWAP" and len(args) == 2:
            (k1, k2) = args
            return _cnot_funcop(nwires, k1, k2) \
                    * _cnot_funcop(nwires, k2, k1) \
                    * _cnot_funcop(nwires, k1, k2)
        elif gate == "QFT" and len(args) == 2:
            (k1, k2) = args
            acc = FuncOperator.identity(is_ketop=True)
            for X in qft_elements(k1, k2):
                acc = Circuit.get_funcoperator(X[0], nwires, X[1:]) * acc
            return acc
        elif gate == "IQFT" and len(args) == 2:
            (k1, k2) = args
            acc = FuncOperator.identity(is_ketop=True)
            for X in iqft_elements(k1, k2):
                acc = Circuit.get_funcoperator(X[0], nwires, X[1:]) * acc
            return acc
        elif gate == "REVERSE" and len(args) == 2:
            (k1, k2) = args
            acc = FuncOperator.identity(is_ketop=True)
            for X in reverse_elements(k1, k2):
                acc = Circuit.get_funcoperator(X[0], nwires, X[1:]) * acc
            return acc
        elif gate == "FUNC" and len(args) == 5 and args[2] == "xyModN":
            (k, q, _, x, N) = args
            return _xymodN_funcop(nwires, k, q, x, N)
        elif gate == "CFUNC" and len(args) == 6 and args[3] == "xyModN":
            (k1, k2, q, _, x, N) = args
            return _cxymodN_funcop(nwires, k1, k2, q, x, N)
        else:
            raise Exception("Invalid circuit element")

    def generate(self, kind: str="mat", **kwargs):
        """
        Generate the single unitary `*Operator` for the circuit.
        """
        if kind == "vec":
            acc = VecOperator.identity(is_ketop=True)
            for X in self.elements:
                acc = Circuit.get_vecoperator(X[0], self.nwires, X[1:]) * acc
            return acc
        elif kind == "mat":
            sparse = kwargs.get("sparse", False)
            acc = MatOperator.identity(self.basis, sparse).action
            for X in self.elements:
                acc = Circuit.get_mat(X[0], self.nwires, X[1:], sparse) @ acc
            return MatOperator(
                acc,
                self.basis,
                is_ketop=True
            )
        elif kind == "func":
            acc = FuncOperator.identity(is_ketop=True)
            for X in self.elements:
                acc = Circuit.get_funcoperator(X[0], self.nwires, X[1:]) * acc
            return acc
        else:
            raise Exception(f"Invalid operator kind: {kind}")

    @staticmethod
    def invert(elements: list[tuple]):
        """
        Reverse a list of quantum circuit elements and replace each constituent
        with its inverse operation.
        """
        _elements = list()
        for X in reversed(elements):
            if X[0] in {"H", "HADAMARD", "CNOT", "NOT", "SWAP", "REVERSE"}:
                _elements.append(X)
            elif X[0] in {"P", "PHASE", "RZ"}:
                (G, k, theta) = X
                _elements.append((G, k, -theta))
            elif X[0] in {"CP", "CPHASE", "CRZ"}:
                (G, k1, k2, theta) = X
                _elements.append((G, k1, k2, -theta))
            elif X[0] == "QFT":
                (G, k1, k2) = X
                _elements.append(("IQFT", k1, k2))
            elif X[0] == "IQFT":
                (G, k1, k2) = X
                _elements.append(("QFT", k1, k2))
            elif X[0] in {"FUNC", "CFUNC"}:
                raise Exception("Inversion of FUNC and CFUNC gates is not implemented")
        return _elements

    def inverted(self):
        """
        Invert the action of self.
        """
        return Circuit(
            self.nwires,
            Circuit.invert(self.elements),
            self.initial_state,
            self.measure
        )

    @staticmethod
    def parse_file(infilename: str):
        """
        Read a list of circuit elements from a given circuit file.
        """
        source = pathlib.Path(infilename)
        with source.open('r') as infile:
            lines = infile.readlines()
        nwires = None
        elements = list()
        measure = False
        initial_state = None
        state_pattern = re.compile(r'\|([01]+)\>')
        gates = {
            "H":        (int,),
            "HADAMARD": (int,),
            "P":        (int, float),
            "PHASE":    (int, float),
            "CNOT":     (int, int),
            "NOT":      (int,),
            "RZ":       (int, float),
            "CRZ":      (int, int, float),
            "CPHASE":   (int, int, float),
            "CP":       (int, int, float),
            "SWAP":     (int, int),
            "QFT":      (int, int),
            "IQFT":     (int, int),
            "REVERSE":  (int, int),
            "FUNC":     (int, int, str, int, int),
            "CFUNC":    (int, int, int, str, int, int),
        }
        for i, line in enumerate(lines):
            largs = line.replace("\n", "").split("#")[0].split(" ")
            try:
                if len(largs) == 0:
                    continue
                elif len(largs) == 1:
                    if largs[0] == "MEASURE":
                        measure = True
                    else:
                        nwires = int(largs[0])
                        initial_state = (
                            StateVec.from_primitives((1.0 + 0.0j, nwires * "0"))
                                .to_statearr(qubit_basis(nwires))
                        )
                elif largs[0] == "INITSTATE" and len(largs) >= 3:
                    if largs[1] == "BASIS" \
                            and (Q := state_pattern.match(largs[2])):
                        nwires = len(Q.group(1))
                        initial_state = (
                            StateVec.from_primitives((1.0 + 0.0j, Q.group(1)))
                                .to_statearr(qubit_basis(nwires))
                        )
                    elif largs[1] == "FILE":
                        F = pathlib.Path(largs[2])
                        try:
                            nwires, initial_state \
                                    = Circuit.parse_initial_state(F)
                        except FileNotFoundError:
                            nwires, initial_state \
                                    = Circuit.parse_initial_state(
                                        source.parent.joinpath(F.name)
                                    )
                        else:
                            raise Exception("Invalid INITSTATE source")
                    else:
                        raise Exception("Invalid INITSTATE source")
                elif largs[0] == "WIRES" and len(largs) >= 2:
                    nwires = int(largs[1])
                    initial_state = (
                        StateVec.from_primitives((1.0 + 0.0j, nwires * "0"))
                            .to_statearr(qubit_basis(nwires))
                    )
                elif largs[0] in gates.keys() \
                        and len(largs) >= len(gates[largs[0]])+1:
                    argtypes = gates[largs[0]]
                    elements.append((
                        largs[0],
                        *[argtype(x) for argtype, x in zip(argtypes, largs[1:])]
                    ))
                else:
                    raise Exception(
                        f"Mal-formed circuit declaration:\n  {line}"
                    )
            except Exception as ERR:
                errmsg = f"Exception occurred while parsing {infilename}:\n" \
                        + "\n".join(e for e in ERR.args)
                raise Exception(errmsg)
        if nwires is None:
            raise Exception("Wire count is undeclared")
        if initial_state is None:
            raise Exception("initial state is undeclared")
        return nwires, elements, initial_state, measure

    @staticmethod
    def parse_initial_state(infilename: str):
        """
        Generate an initial state from a given initial state file.
        """
        source = pathlib.Path(infilename)
        with source.open('r') as infile:
            lines = [line.split(" ") for line in infile.readlines()]
        N = len(lines)
        if np.log2(N) - np.floor(np.log2(N)) > 0:
            raise Exception("Number of components must be a power of 2")
        nwires = int(np.log2(N))
        try:
            initial_state = StateArr(
                np.array([complex(float(line[0]), float(line[1])) \
                        for line in lines]),
                qubit_basis(nwires)
            ).normalized()
        except ValueError:
            raise Exception("Mal-formed initial state component")
        return nwires, initial_state

    def to_atomic(self):
        """
        Replace all circuit elements with equivalent sub-routines using only
        H/P/CNOT gates.
        """
        _elements = list()
        for X in self.elements:
            if X[0] in {"H", "HADAMARD"}:
                _elements.append(X)
            elif X[0] == {"P", "PHASE"}:
                _elements.append(X)
            elif X[0] == "CNOT":
                _elements.append(X)
            elif X[0] == "NOT":
                (_, k) = X
                _elements.append(("H", k))
                _elements.append(("P", k, np.pi))
                _elements.append(("H", k))
            elif X[0] == "RZ":
                (_, k, theta) = X
                _elements.append(("P", k, theta/2))
                _elements.append(("H", k))
                _elements.append(("P", k, np.pi))
                _elements.append(("H", k))
                _elements.append(("P", k, -theta/2))
                _elements.append(("H", k))
                _elements.append(("P", k, np.pi))
                _elements.append(("H", k))
            elif X[0] == "CRZ":
                (_, k1, k2, theta) = X
                _elements.append(("P", k2, theta/2))
                _elements.append(("CNOT", k1, k2))
                _elements.append(("P", k2, -theta/2))
                _elements.append(("CNOT", k1, k2))
            elif X[0] in {"CP", "CPHASE"}:
                (_, k1, k2, theta) = X
                _elements.append(("P", k2, theta/2))
                _elements.append(("CNOT", k1, k2))
                _elements.append(("P", k2, -theta/2))
                _elements.append(("CNOT", k1, k2))
                _elements.append(("P", k1, theta/2))
            elif X[0] == "SWAP":
                (_, k1, k2) = X
                _elements.append(("CNOT", k1, k2))
                _elements.append(("CNOT", k2, k1))
                _elements.append(("CNOT", k1, k2))
            elif X[0] == "QFT":
                (_, k1, k2) = X
                _elements += qft_elements(k1, k2)
            elif X[0] == "IQFT":
                (_, k1, k2) = X
                _elements += iqft_elements(k1, k2)
            elif X[0] == "REVERSE":
                (_, k1, k2) = X
                _elements += reverse_elements(k1, k2)
            elif X[0] == "FUNC":
                _elements.append(X)
            elif X[0] == "CFUNC":
                _elements.append(X)
        return Circuit(self.nwires, _elements, self.initial_state, self.measure)

    @staticmethod
    def from_file(infilename: str):
        """
        Generate a `Circuit` from a given circuit file.
        """
        return Circuit(*Circuit.parse_file(infilename))

    def to_file(self, outfilename: str):
        """
        Write a description of self to a file.
        """
        target = pathlib.Path(outfilename)
        outfile = target.open('w')
        outfile.write(f"{self.nwires}\n")
        if len(S := self.initial_state.to_statevec()) > 1:
            starget = target.parent.joinpath(target.stem + "_initstate.txt")
            with starget.open('w') as sfile:
                for b, a in self.initial_state:
                    sfile.write(f"{a.real:+.8f} {a.imag:+.8f}\n")
            outfile.write(f"INITSTATE FILE {starget}\n")
        else:
            s = list(S.components.keys())[0].label
            outfile.write(f"INITSTATE BASIS |{s}>\n")
        for X in self.elements:
            outfile.write(" ".join(str(x) for x in X) + "\n")
        if self.measure:
            outfile.write("MEASURE")
        outfile.close()

    @staticmethod
    def random(Ngates: int, nwires: int, measure: bool=False):
        """
        Generate a `Circuit` with randomly chosen operations.
        """
        gates = ["H", "P", "CNOT"]
        elements = list()
        for k in range(Ngates):
            X = random.choice(gates)
            if X == "H":
                elements.append((X, random.randrange(0, nwires)))
            elif X == "P":
                elements.append(
                    (X, random.randrange(0, nwires), 2*np.pi*random.random())
                )
            elif X == "CNOT":
                k1 = random.randrange(0, nwires)
                k2 = k1
                while k2 == k1:
                    k2 = random.randrange(0, nwires)
                elements.append((X, k1, k2))
        initial_state = (
            StateVec.from_primitives((1.0 + 0.0j, nwires * "0"))
                .to_statearr(qubit_basis(nwires))
        )
        return Circuit(nwires, elements, initial_state, measure)

    def run(self, kind: str="mat", stepwise: bool=False, **kwargs):
        """
        Run self to calculate an output state and return a `CircuitOutput`.
        """
        sparse = kwargs.get("sparse", False)
        if kind == "vec":
            state = self.initial_state.to_statevec()
            if stepwise:
                for X in self.elements:
                    state = Circuit.get_vecoperator(
                        X[0],
                        self.nwires,
                        X[1:]
                    ) * state
            else:
                state = self.generate(kind, **kwargs) * state
        elif kind == "mat":
            state = self.initial_state.to_statearr(self.basis)
            if stepwise:
                for X in self.elements:
                    state = Circuit.get_matoperator(
                        X[0],
                        self.nwires,
                        X[1:],
                        sparse
                    ) * state
            else:
                state = self.generate(kind, **kwargs) * state
        elif kind == "func":
            state = self.initial_state.to_statevec()
            if stepwise:
                for X in self.elements:
                    state = Circuit.get_funcoperator(
                        X[0],
                        self.nwires,
                        X[1:]
                    ) * state
            else:
                state = self.generate(kind, **kwargs) * state
        else:
            raise Exception(f"Invalid operator kind: {kind}")
        if isinstance(self.measure, bool) and self.measure:
            return CircuitOutput(
                state,
                self.basis,
                state.measure(2**self.nwires * 1000, qubit_to_int)
            )
        elif isinstance(self.measure, int):
            return Circuitoutput(
                state,
                self.basis,
                state.measure(self.measure, qubit_to_int)
            )
        else:
            return CircuitOutput(
                state,
                self.basis,
                None
            )

