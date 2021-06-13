from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as pp
import matplotlib.widgets as wig
import lib.plotdefs as pd
import random
import copy

"""
Provides basic machinery for working with quantum states and operators. It is
assumed that we only work with a single basis of orthonormal states, and that
outer products are disallowed (multiplication of states can give only scalar
values or errors).

States are represented either as `StateVec`s, which are dicts mapping basis
states to amplitudes, or `StateArr`s, which are numpy arrays containing the
full decomposition of the state in a particular basis. Either form can be
converted to the other.

Operators are represented (roughly in order of computing speed) as
`MatOperator`s, `VecOperator`s, or `FuncOperator`s. `MatOperator`s contain the
full matrix representation of an operator and only multiply `StateArr`s;
`VecOperator`s map basis states to an appropriate output `StateVec` and only
multiply `StateVec`s. Either of these forms can be converted to the other.
`FuncOperator`s use pre-defined Python functions to perform actions (and hence
are the only operators here that can operate on an infinite basis) and only
multiply `StateVec`s, but can be converted to either of the other two
representations.

`MatOperator`s and `StateArr`s can be stored as sparse arrays using scipy.sparse
to save memory.
"""

Number = (int, float, complex)
Function = type(lambda:0)

class BasisState:
    """
    Represents the most basic form of a basis state (i.e. neither a ket nor a
    bra; really just a label).
    """
    def __init__(self, label: str):
        self.label = label

    def __eq__(self, other: BasisState):
        assert isinstance(other, BasisState)
        return self.label == other.label

    def __hash__(self):
        return self.label.__hash__()

    def __str__(self):
        return f"`{self.label}`"

    def __repr__(self):
        return f"BasisState(`{self.label}`)"

    def __gt__(self, other):
        if isinstance(other, BasisState):
            return self.label > other.label
        else:
            raise NotImplemented

    def __ge__(self, other):
        if isinstance(other, BasisState):
            return self.label >= other.label
        else:
            raise NotImplemented

    def __lt__(self, other):
        if isinstance(other, BasisState):
            return self.label < other.label
        else:
            raise NotImplemented

    def __le__(self, other):
        if isinstance(other, BasisState):
            return self.label <= other.label
        else:
            raise NotImplemented

    def __eq__(self, other):
        if isinstance(other, BasisState):
            return self.label == other.label
        else:
            raise NotImplemented

    def __neq__(self, other):
        if isinstance(other, BasisState):
            return self.label != other.label
        else:
            raise NotImplemented

    def to_statevec(self, is_ket: bool=True):
        return StateVec({self: 1.0 + 0.0j}, is_ket)

    def to_statearr(self, basis: Basis, is_ket: bool=True, sparse: bool=False):
        if sparse:
            if is_ket:
                _components = sp.csr_matrix(
                    [[(1.0 + 0.0j) if self == b else (0.0 + 0.0j)]
                        for b in basis]
                )
            else:
                _components = sp.csr_matrix(
                    [[(1.0 + 0.0j) if self == b else (0.0 + 0.0j)
                        for b in basis]]
                )
        else:
            _components = np.array(
                [(1.0 + 0.0j) if self == b else (0.0 + 0.0j)
                    for b in basis]
            )
        return StateArr(_components, basis, is_ket)

    def relabeled(self, label_func):
        return BasisState(label_func(self.label))

class Basis:
    """
    List of `BasisState`s.
    """
    def __init__(self, *states: tuple[BasisState]):
        self.states = states

    @staticmethod
    def from_labels(*labels: tuple[str]):
        return Basis(*[BasisState(l) for l in labels])

    def __eq__(self, other: Basis):
        assert isinstance(other, Basis)
        return self.states == other.states

    def __iter__(self):
        return iter(self.states)

    def __getitem__(self, pos: int):
        return self.states[pos]

    def __contains__(self, elem: BasisState):
        assert isinstance(elem, BasisState)
        return elem in self.states

    def index(self, value):
        return self.states.index(value)

    def __str__(self):
        return "{"+", ".join(str(s) for s in self.states)+"}"

    def __repr__(self):
        return "Basis("+", ".join(str(s) for s in self.states)+")"

    def __len__(self):
        return len(self.states)

    def relabeled(self, label_func):
        return Basis(*[s.relabeled(label_func) for s in self.states])

class StateVec:
    """
    Quantum state represented as a collection of Dirac kets (or bras). The state
    is internally stored as dict[BasisState: complex] (and is indexed as such).
    Interacts only with other `StateVec`s, `VecOperator`s, and `FuncOperator`s.
    """
    def __init__(self, components: dict[BasisState, complex],
            is_ket: bool=True):
        self.components = {b: complex(a) for b, a in components.items()}
        self.is_ket = is_ket

    @staticmethod
    def from_states(*components: tuple[StateVec]):
        """
        Construct a single `StateVec` from many other `StateVec`s, summing
        amplitudes for identical basis states and removing any zero-amplitude
        terms in the process. All inputs to this function must have the same
        `.is_ket` value.
        """
        kind_sum = sum(map(lambda S: S.is_ket, components))
        assert kind_sum in {0, len(components)}
        _is_ket = bool(kind_sum)
        _components = dict()
        for S in components:
            for b, a in S:
                if a != 0:
                    _components[b] = _components[b] + complex(a) \
                            if b in _components.keys() else complex(a)
                    if abs(_components[b]) == 0:
                        del _components[b]
        return StateVec(_components, _is_ket)

    @staticmethod
    def from_primitives(*components: tuple[complex, str], is_ket: bool=True):
        """
        Construct a single `StateVec` from many tuple[complex, str] pairs. Does
        not combine like terms or remove zero-amplitude terms.
        """
        return StateVec(
            {BasisState(l): a for a, l in components if a != 0},
            is_ket
        )

    @staticmethod
    def zero(is_ket: bool=True):
        return StateVec(dict(), is_ket)

    def conjugate(self):
        return StateVec(
            {b: a.conjugate() for b, a in self},
            not self.is_ket
        )

    def hc(self):
        return self.conjugate()

    def to_statevec(self):
        """
        Deep copy of self.
        """
        return copy.deepcopy(self)

    def to_statearr(self, basis: Basis, sparse=False):
        """
        Return a `StateArr` containing the decomposition of self in the given
        `Basis`.
        """
        if sparse:
            if self.is_ket:
                _components = sp.csr_matrix([[self.__mul__(b)] for b in basis])
            else:
                _components = sp.csr_matrix([[self.__mul__(b) for b in basis]])
        else:
            _components = np.array([self.__mul__(b) for b in basis])
        return StateArr(_components, basis, self.is_ket)

    def keys(self):
        return self.components.keys()

    def get(self, key, default):
        return self.components.get(key, default)

    def __getitem__(self, state: BasisState):
        assert isinstance(state, BasisState)
        return self.components.get(state, 0.0 + 0.0j)

    def __contains__(self, key):
        return key in self.components

    def __neg__(self):
        return StateVec({b: a.__neg__() for b, a in self}, self.is_ket)

    def __add__(self, other: StateVec):
        assert isinstance(other, StateVec) and self.is_ket == other.is_ket
        return StateVec.from_states(self, other)

    def __radd__(self, other: StateVec):
        return other.__add__(self)

    def __sub__(self, other: StateVec):
        return self.__add__(other.__neg__())

    def __rsub__(self, other: StateVec):
        return other.__sub__(self)

    def __mul__(self,
            other: BasisState | StateVec | VecOperator | FuncOperator | Number):
        assert isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateVec)
                    and (not self.is_ket and other.is_ket)) \
                or (isinstance(other, VecOperator)
                    and (not self.is_ket and not other.is_ketop)) \
                or (isinstance(other, FuncOperator)
                    and (not self.is_ket and not other.is_ketop))
        if isinstance(other, BasisState):
            return self.components.get(other, 0.0 + 0.0j)
        elif isinstance(other, StateVec):
            shared_keys = set(self.keys()).intersection(set(other.keys()))
            return sum([self[b] * other[b] for b in shared_keys])
        elif isinstance(other, Number):
            return StateVec({b: a * other for b, a in self}, self.is_ket)
        else:
            return other.__rmul__(self)

    def __rmul__(self,
            other: BasisState | StateVec | VecOperator | FuncOperator | Number):
        assert isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateVec)
                    and (self.is_ket and not other.is_ket)) \
                or (isinstance(other, VecOperator)
                    and (self.is_ket and other.is_ketop)) \
                or (isinstance(other, FuncOperator)
                    and (self.is_ket and other.is_ketop))
        if isinstance(other, BasisState):
            return self.components.get(other, 0.0 + 0.0j)
        elif isinstance(other, StateVec):
            shared_keys = set(self.keys()).intersection(set(other.keys()))
            return sum([self[b] * other[b] for b in shared_keys])
        elif isinstance(other, Number):
            return StateVec({b: other * a for b, a in self}, self.is_ket)
        else:
            return other.__mul__(self)

    def __truediv__(self, other: Number):
        assert isinstance(other, Number)
        return StateVec({b: a / other for b, a in self})

    def __str__(self):
        if len(self) == 0:
            return "0"
        else:
            return " + ".join(f"({a:+.5f})" \
                    + ("|" if self.is_ket else "<") \
                    + f"{b.label}" \
                    + (">" if self.is_ket else "|")
                for b, a in self if abs(a) > 1e-12
            )

    def __repr__(self):
        return "StateVec({" \
                + ", ".join(f"{b}: {a:g}" for b, a in self) \
                + "})"

    def __len__(self):
        return len(self.components.keys())

    def __iter__(self):
        return iter(self.components.items())

    def relabeled(self, label_func):
        return StateVec(
            {b.relabeled(label_func): a for b, a in self},
            self.is_ket
        )

    def norm(self):
        return np.sqrt(sum(abs(a)**2 for a in self.components.values()))

    def normalized(self):
        """
        Return a copy of self where amplitudes are properly normalized.
        """
        return self.__truediv__(self.norm())

    def measure_single(self, state_func=None):
        """
        Perform a single measurement on self, returning any constituent
        `BasisState` with appropriate probabilities. A function taking as input
        a single `BasisState` and returning anything can optionally be provided
        to call on the result of the measurement.
        """
        N2 = self.norm()**2
        r = random.random()
        acc = 0
        for s, a in self.components.items():
            acc += abs(a)**2/N2
            if r < acc:
                return state_func(s) if state_func is not None else s

    def measure(self, N: int=1000, state_func=None):
        """
        Perform `measure_single` on self for `N` times and return the results in
        a list.
        """
        return [self.measure_single(state_func) for k in range(N)]

    def most_probable(self):
        """
        Return the most probable result of a single measurement.
        """
        P = 0.0
        B = list(self.keys())[0]
        for b, a in self:
            if abs(a)**2 > P:
                P = abs(a)**2
                B = b
        return B

    def draw(self, levels: list[BasisState]=None, plotter: pd.Plotter=None):
        """
        Draw the Bloch-sphere representation of self. Requires that self be a
        superposition of a two-level system or that `levels` be specified with
        two `BasisState`s.
        """
        if levels is None:
            assert len(self) == 2
            b0, b1 = sorted(self.keys())
            alpha = complex(self[b0])
            beta = complex(self[b1])
        else:
            assert len(levels) == 2
            assert sum(map(lambda b: b in self, levels)) > 0
            b0, b1 = levels
            alpha = complex(self.get(b0, 0.0 + 0.0j))
            beta = complex(self.get(b1, 0.0 + 0.0j))
        l = "|" if self.is_ket else "\\rangle"
        r = "\\rangle" if self.is_ket else "|"
        if plotter is None:
            P = _gen_bloch_sphere([
                "$\\left" + l + b.label + "\\right" + r + "$"
                for b in [b0, b1]
            ])
        else:
            P = plotter
        return _draw_bloch(alpha, beta, self.is_ket, P)

class StateArr:
    """
    Quantum state represented as a collection of complex amplitudes giving a
    decomposition in a particular `Basis`. The state is internally stored as
    numpy.ndarray[complex128] or scipy.sparse.csr_matrix[complex128] (and is
    indexed as such). Interacts only with other `StateArr`s or `MatOperator`s.
    """
    def __init__(self, components: np.ndarray | sp.csr_matrix, basis: Basis,
            is_ket: bool=True):
        assert components.shape == (len(basis),) \
                or components.shape == (len(basis), 1) \
                or components.shape == (1, len(basis))
        self.components = components.astype(complex)
        self.basis = basis
        self.is_ket = is_ket
        self.is_sparse = not isinstance(components, np.ndarray)

    @staticmethod
    def zero(basis: Basis, is_ket: bool=True, sparse: bool=False):
        if sparse:
            _components = sp.csr_matrix(
                len(basis) * [[0]] if is_ket else len(basis) * [0]
            )
        else:
            _components = np.zeros(len(basis))
        return StateArr(_components, basis, is_ket)

    def conjugate(self):
        return StateArr(
            self.components.conjugate().T,
            self.basis,
            not self.is_ket
        )

    def hc(self):
        return self.conjugate()

    def to_statevec(self):
        """
        Return a `StateVec` containing the minimal equivalent ket-bra
        representation of self (i.e. with zero components removed).
        """
        return StateVec({b: a for b, a in self}, self.is_ket)

    def to_statearr(self, basis: Basis):
        """
        Deep copy of self.
        """
        return copy.deepcopy(self)

    def to_sparse(self):
        """
        Convert to a sparse internal representation.
        """
        if self.is_sparse:
            return self
        else:
            return StateArr(
                sp.csr_matrix(self.components),
                self.basis,
                self.is_ket
            )

    def to_dense(self):
        """
        Convert to a dense internal representation.
        """
        if self.is_sparse:
            return StateArr(
                self.components.toarray(),
                self.basis,
                self.is_ket
            )
        else:
            return self

    def __getitem__(self, idx):
        return self.components[idx]

    def __contains__(self, b):
        return b in self.basis

    def __neg__(self):
        return StateArr(self.components.__neg__(), self.basis, self.is_ket)

    def __add__(self, other: StateArr):
        assert isinstance(other, StateArr)
        assert self.is_ket == other.is_ket
        assert self.basis == other.basis
        return StateArr(
            self.components.__add__(other.components),
            self.basis,
            self.is_ket
        )

    def __radd__(self, other: StateArr):
        return self.__add__(other)

    def __sub__(self, other: StateArr):
        return self.__add__(other.__neg__())

    def __rsub__(self, other: StateArr):
        return other.__sub__(self)

    def __mul__(self, other: BasisState | StateArr | MatOperator | Number):
        assert isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateArr)
                    and not self.is_ket and other.is_ket
                    and self.basis == other.basis) \
                or (isinstance(other, MatOperator)
                    and not self.is_ket and not other.is_ketop
                    and self.basis == other.basis)
        if isinstance(other, BasisState):
            return (0.0 + 0.0j) if other not in self.basis \
                    else self.components[self.basis.index(other)]
        elif isinstance(other, StateArr):
            return self.components.__mul__(other.components).sum()
        elif isinstance(other, Number):
            return StateArr(
                self.components.__mul__(other),
                self.basis,
                self.is_ket
            )
        else:
            return other.__rmul__(self)

    def __rmul__(self, other: BasisState | StateArr | MatOperator | Number):
        assert isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateArr)
                    and self.is_ket and not other.is_ket
                    and self.basis == other.basis) \
                or (isinstance(other, MatOperator)
                    and self.is_ket and other.is_ketop
                    and self.basis == other.basis)
        if isinstance(other, BasisState):
            return (0.0 + 0.0j) if other not in self.basis \
                    else self.components[self.basis.index(other)]
        elif isinstance(other, StateArr):
            return self.components.__mul__(other.components).sum()
        elif isinstance(other, Number):
            return StateArr(
                self.components.__mul__(other),
                self.basis,
                self.is_ket
            )
        else:
            return other.__mul__(self)

    def __truediv__(self, other: Number):
        assert isinstance(other, Number)
        return StateArr(
            self.components.__truediv__(other),
            self.basis,
            self.is_ket
        )

    def __str__(self):
        return "\n".join(f"[ {a:+.5f} ] " \
                + ("|" if self.is_ket else "<") \
                + f"{b.label}" \
                + (">" if self.is_ket else "|")
            for b, a in self
        )

    def __repr__(self):
        return "StateArr(\n" \
                + "  "+"["+", ".join(f"{a:g}" for a in self.components)+"],\n" \
                + "  "+self.basis.__repr__()+",\n" \
                + "  "+f"is_ket={self.is_ket}\n" \
                + ")"

    def __len__(self):
        return len(self.components)

    def __iter__(self):
        return iter(zip(self.basis, self.components))

    def relabeled(self, label_func):
        return StateArr(
            self.components,
            self.basis.relabeled(label_func),
            self.is_ket
        )

    def norm(self):
        return np.sqrt(sum(abs(a)**2 for a in self.components))

    def normalized(self):
        """
        Return a copy of self where amplitudes are properly normalized.
        """
        return self.__truediv__(self.norm())

    def measure_single(self, state_func=None):
        """
        Perform a single measurement on self, returning any constituent
        `BasisState` with appropriate probabilities. A function taking as input
        a single `BasisState` and returning anything can optionally be provided
        to call on the result of the measurement.
        """
        N2 = self.norm()**2
        r = random.random()
        acc = 0
        for s, a in zip(self.basis, self.components):
            acc += abs(a)**2/N2
            if r < acc:
                return state_func(s) if state_func is not None else s

    def measure(self, N: int=1000, state_func=None):
        """
        Perform `measure_single` on self for `N` times and return the results in
        a list.
        """
        return [self.measure_single(state_func) for k in range(N)]

    def most_probable(self):
        """
        Return the most probable result of a single measurement.
        """
        return self.basis[np.argmax(np.abs(self.components)**2)]

    def draw(self, levels: list[BasisState]=None, plotter: pd.Plotter=None):
        """
        Draw the Bloch-sphere representation of self. Requires that self be a
        superposition of a two-level system or that `levels` be specified with
        two `BasisState`s.
        """
        if levels is None:
            assert len(self) == 2
            b0, b1 = self.basis
            alpha, beta = self.components[0], self.components[1]
        else:
            assert len(levels) == 2
            assert sum(map(lambda b: b in self.basis, levels)) > 0
            b0, b1 = levels
            alpha = self.components[self.basis.index(b0)] \
                    if b0 in self.basis else (0.0 + 0.0j)
            beta = self.components[self.basis.index(b1)] \
                    if b1 in self.basis else (0.0 + 0.0j)
        l = "|" if self.is_ket else "\\rangle"
        r = "\\rangle" if self.is_ket else "|"
        if plotter is None:
            P = _gen_bloch_sphere([
                "$\\left" + l + b.label + "\\right" + r + "$"
                for b in [b0, b1]
            ])
        else:
            P = plotter
        return _draw_bloch(alpha, beta, self.is_ket, P)

class VecOperator:
    """
    Quantum operator represented as a collection of `BasisState`-to-`StateVec`
    mappings. The action of the operator is internally stored as
    dict[BasisState: StateVec] (and is indexed as such). Interacts only with
    other `VecOperator`s, `FuncOperator`s, and `StateVec`s.
    """
    def __init__(self, action: dict[BasisState, StateVec], is_ketop: bool=True,
            default_zero: bool=True):
        self.action = action
        self.is_ketop = is_ketop
        self.default_zero = default_zero

    @staticmethod
    def identity(is_ketop=True):
        return VecOperator(
            dict(),
            is_ketop,
            default_zero=False,
        )

    def adjoint(self):
        return VecOperator(
            {b: self.action[b].conjugate() for b in self.keys()},
            not self.is_ketop,
            self.default_zero,
        )

    def hc(self):
        return self.adjoint()

    def to_matoperator(self, basis: Basis, sparse: bool=False):
        """
        Convert to an equivalent `MatOperator` representation in the given
        `Basis`.
        """
        return MatOperator(
            (sp.csr_matrix if sparse else np.array)([
                [self[B].__rmul__(b) if self.is_ketop
                    else self[B].__mul__(b).conjugate() for b in basis]
                for B in basis
            ]).T,
            basis,
            self.is_ketop
        )

    def keys(self):
        return self.action.keys()

    def get(self, key, default):
        return self.action.get(key, default)

    def __getitem__(self, b: BasisState):
        assert isinstance(b, BasisState)
        return self.action.get(
            b,
            StateVec.zero(self.is_ketop) if self.default_zero \
                    else StateVec({b: 1.0 + 0.0j}, self.is_ketop)
        )

    def __neg__(self):
        return VecOperator(
            {b: -s for b, s in self.action.keys()},
            self.is_ketop,
            self.default_zero,
        )

    def __add__(self, other: VecOperator | FuncOperator | Number):
        assert (isinstance(other, VecOperator)
                    and self.is_ketop == other.is_ketop) \
                or (isinstance(other, FuncOperator)
                    and self.is_ketop == other.is_ketop) \
                or isinstance(other, Number)
        if isinstance(other, VecOperator):
            all_keys = set(self.keys()).union(set(other.keys()))
            return VecOperator(
                {b: self[b] + other[b] for b in all_keys},
                self.is_ketop,
                default_zero=True,
            )
        elif isinstance(other, Number):
            return self.__add__(VecOperator.identity().__mul__(other))
        else:
            return other.__radd__(self)

    def __radd__(self, other: VecOperator | FuncOperator | Number):
        assert (isinstance(other, VecOperator)
                    and self.is_ketop == other.is_ketop) \
                or (isinstance(other, FuncOperator)
                    and self.is_ketop == other.is_ketop) \
                or isinstance(other, Number)
        if isinstance(other, VecOperator):
            return other.__add__(self)
        elif isinstance(other, Number):
            return VecOperator.identity().__mul__(other).__add__(self)
        else:
            return other.__add__(self)

    def __sub__(self, other: VecOperator | FuncOperator | Number):
        if isinstance(other, VecOperator):
            return self.__add__(other.__neg__())
        elif isinstance(other, Number):
            return self.__add__(VecOperator.identity().__mul__(-other))
        else:
            return other.__neg__().__radd__(self)

    def __rsub__(self, other: VecOperator | FuncOperator | Number):
        if isinstance(other, VecOperator):
            return other.__add__(self.__neg__())
        elif isinstance(other, Number):
            return VecOperator.identity().__mul__(other).__add__(self.__neg__())
        else:
            return other.__add__(self.__neg__())

    def __mul__(self,
            other: BasisState | StateVec | VecOperator | FuncOperator | Number):
        assert isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateVec)
                    and self.is_ketop and other.is_ket) \
                or (isinstance(other, (VecOperator, FuncOperator))
                    and self.is_ketop == other.is_ketop)
        if isinstance(other, BasisState):
            return self[other]
        elif isinstance(other, StateVec):
            return StateVec.from_states(*[
                self[b].__rmul__(a) for b, a in other
            ])
        elif isinstance(other, VecOperator):
            all_keys = set(self.keys()).union(set(other.keys()))
            return VecOperator(
                {b: self.__mul__(other[b]) for b in all_keys},
                self.is_ketop,
                self.default_zero and other.default_zero,
            )
        elif isinstance(other, Number):
            return VecOperator(
                {b: s__mul__(other) for b, s in self.action.items()},
                self.is_ketop,
                self.default_zero,
            )
        else:
            return other.__rmul__(self)

    def __rmul__(self,
            other: BasisState | StateVec | VecOperator | FuncOperator | Number):
        assert (isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateVec)
                    and not self.is_ketop and not other.is_ket) \
                or (isinstance(other, (VecOperator, FuncOperator))
                    and self.is_ketop == other.is_ketop))
        if isinstance(other, BasisState):
            return self[other]
        elif isinstance(other, StateVec):
            return StateVec.from_states(*[
                self[b].__rmul__(a) for b, a in other
            ])
        elif isinstance(other, VecOperator):
            all_keys = set(self.keys()).union(set(other.keys()))
            return VecOperator(
                {b: self.__rmul__(other[b]) for b in all_keys},
                self.is_ketop,
                other.default_zero and self.default_zero,
            )
        elif isinstance(other, Number):
            return VecOperator(
                {b: s.__rmul__(other) for b, s in self.action.keys()},
                self.is_ketop,
                self.default_zero,
            )
        else:
            return other.__mul__(self)

    def __truediv__(self, other: Number):
        assert isinstance(other, Number)
        return VecOperator(
            {b: s.__truediv__(other) for b, s in self.action.items()},
            self.is_ketop,
            self.default_zero,
        )

    def __iter__(self):
        return iter(self.action.items())

    def __str__(self):
        return "VecOperator({\n" \
                + (
                    "\n".join(f"{b} -> {s}"
                    for b, s in self.action.items())
                ) \
                + f"\n}}, default_zero = {self.default_zero})"

    def relabeled(self, label_func):
        return VecOperator(
            {b.relabeled(label_func): s.relabeled(label_func) for b, s in self},
            self.is_ketop,
            self.default_zero,
        )

class MatOperator:
    """
    Quantum operator represented as a an array of matrix elements for a given
    `Basis`. The action of the operator is internally stored as
    numpy.ndarray[complex128] or scipy.sparse.csr_matrix[complex128] (and is
    indexed as such). Interacts only with other `MatOperator`s and `StateVec`s.
    """
    def __init__(self, action: np.ndarray | sp.csr_matrix, basis: Basis,
            is_ketop: bool=True):
        assert action.shape == (len(basis), len(basis))
        self.action = action
        self.basis = basis
        self.is_ketop = is_ketop
        self.is_sparse = False if isinstance(action, np.ndarray) else True

    @staticmethod
    def identity(basis: Basis, is_ketop: bool=True, sparse: bool=False):
        _action = sp.identity(len(basis), dtype=complex, format="csr") \
                if sparse else np.identity(len(basis), dtype=complex)
        return MatOperator(_action, basis, is_ketop)

    def adjoint(self):
        return MatOperator(
            self.action.conjugate().T,
            self.basis,
            not self.is_ketop
        )

    def hc(self):
        return self.adjoint()

    def to_vecoperator(self):
        """
        Convert to an equivalent `VecOperator` representation.
        """
        _action = dict()
        for i in range(len(self.basis)):
            _action[self.basis[i]] = StateArr(
                self.action[:,i] if self.is_ketop else self.action[i,:],
                self.basis,
                self.is_ketop
            ).to_statevec()
        return VecOperator(
            _action,
            self.is_ketop,
            default_zero=True,
        )

    def to_sparse(self):
        """
        Convert to a sparse internal representation.
        """
        if self.is_sparse:
            return self
        else:
            return MatOperator(
                sp.csr_matrix(self.action),
                self.basis,
                self.is_ketop
            )

    def to_dense(self):
        """
        Convert to a dense internal representation.
        """
        if self.is_sparse:
            return MatOperator(
                self.action.toarray(),
                self.basis,
                self.is_ketop
            )
        else:
            return self

    def __neg__(self):
        return MatOperator(
            self.action__neg__(),
            self.basis,
            self.is_ketop
        )

    def __add__(self, other: MatOperator | Number):
        assert isinstance(other, (MatOperator, *Number))
        if isinstance(other, MatOperator) and self.basis == other.basis:
            return MatOperator(
                self.action.__add__(other.action),
                self.basis,
                self.is_ketop
            )
        elif isinstance(other, Number):
            return MatOperator(
                self.action.__add__(other),
                self.basis,
                self.is_ketop
            )
        else:
            raise Exception("Invalid add")

    def __sub__(self, other: MatOperator | Number):
        return self.__add__(other.__neg__())

    def __mul__(self, operand: StateArr | MatOperator | Number):
        assert isinstance(operand, (StateArr, MatOperator, *Number))
        if isinstance(operand, StateArr) \
                and operand.is_ket and self.basis == operand.basis:
            return StateArr(
                self.action.__matmul__(operand.components),
                self.basis,
                operand.is_ket
            )
        elif isinstance(operand, MatOperator) and self.basis == operand.basis:
            return MatOperator(
                self.action.__matmul__(operand.action),
                self.basis,
                self.is_ketop
            )
        elif isinstance(operand, Number):
            return MatOperator(
                self.action.__mul__(operand),
                self.basis,
                self.is_ketop
            )
        else:
            raise Exception("Invalid multiply")

    def __rmul__(self, operand: StateArr | MatOperator | Number):
        assert isinstance(operand, (StateArr, MatOperator, *Number))
        if isinstance(operand, StateArr) \
                and not operand.is_ket and self.basis == operand.basis:
            return StateArr(
                operand.components.__matmul__(self.action),
                self.basis,
                operand.is_ket
            )
        elif isinstance(operand, MatOperator) and self.basis == operand.basis:
            return MatOperator(
                operand.action.__matmul__(self.action),
                self.basis,
                other.is_ketop
            )
        elif isinstance(operand, Number):
            return MatOperator(
                self.action.__mul__(operand),
                self.basis,
                self.is_ketop
            )
        else:
            raise Exception("Invalid multiply")

    def __truediv__(self, other: Number):
        assert isinstance(other, Number)
        return MatOperator(
            self.action.__truediv__(other),
            self.basis,
            self.is_ketop
        )

    def __str__(self):
        return "MatOperator(\n" \
                + str(self.action) + "\n" \
                + f") basis = {str(self.basis)}"

    def relabeled(self, label_func):
        return MatOperator(
            self.action,
            self.basis.relabeled(label_func),
            self.is_ketop
        )

class FuncOperator:
    """
    Quantum operator represented as a black bax mapping `BasisState`s to
    corresponding `StateVec`s. The action of the operator is internally stored
    as a Python function. Interacts only with all `*Operator`s and `StateVec`s.

    Functions must have the signature
        action(BasisState) -> StateVec(is_ket=is_ketop)
    upon construction using __init__. Combinations of `FuncOperator`s with other
    `*Operator`s are created by composing function calls (and always return new
    `FuncOperator`s), so be wary of having too many such interactions.
    """
    def __init__(self, action: Function, is_ketop=True):
        self.action = action
        self.is_ketop = is_ketop

    @staticmethod
    def identity(is_ketop=True):
        return FuncOperator(
            lambda b: StateVec({b: 1.0 + 0.0j}, is_ketop),
            is_ketop
        )

    def adjoint(self):
        return FuncOperator(
            lambda b: self.action(b).adjoint(),
            not self.is_ketop
        )

    def hc(self):
        return self.adjoint()

    def to_matoperator(self, basis: Basis, sparse: bool=False):
        """
        Convert to an equivalent `MatOperator` representation in the given
        `Basis`.
        """
        return MatOperator(
            (sp.csr_matrix if sparse else np.array)([
                [(self * B).__rmul__(b) if self.is_ketop
                    else (B * self).__mul__(b).conjugate() for b in basis]
                for B in basis
            ]).T,
            basis,
            self.is_ketop
        )

    def to_vecoperator(self, basis: Basis, default_zero: bool=True):
        """
        Convert to an equivalent `VecOperator` representation computed for a
        given `Basis`.
        """
        return VecOperator(
            {b: self.action(b) for b in basis},
            self.is_ketop,
            default_zero,
        )

    def __neg__(self):
        return FuncOperator(
            lambda b: self.action(b).__neg__(),
            self.is_ketop
        )

    def __add__(self, other: VecOperator | FuncOperator | Number):
        assert (isinstance(other, VecOperator)
                    and self.is_ketop == other.is_ketop) \
                or (isinstance(other, FuncOperator)
                    and self.is_ketop == other.is_ketop) \
                or isinstance(other, Number)
        if isinstance(other, VecOperator):
            return FuncOperator(
                lambda b: self.action(b) \
                        + ((other * b) if other.is_ketop else (b * other)),
                self.is_ketop
            )
        elif isinstance(other, FuncOperator):
            return FuncOperator(
                lambda b: self.action(b) + other.action(b),
                self.is_ketop
            )
        elif isinstance(other, Number):
            return FuncOperator(
                lambda b: self.action(b) + StateVec({b: other}, self.is_ketop),
                self.is_ketop
            )

    def __radd__(self, other: VecOperator | FuncOperator | Number):
        return self.__add__(other)

    def __sub__(self, other: VecOperator):
        return self.__add__(other.__neg__())

    def __rsub__(self, other: VecOperator):
        return self.__neg__().__add__(other)

    def __mul__(self,
            other: BasisState | StateVec | VecOperator | FuncOperator | Number):
        assert isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateVec)
                    and self.is_ketop and other.is_ket) \
                or (isinstance(other, (VecOperator, FuncOperator))
                    and self.is_ketop == other.is_ketop)
        if isinstance(other, BasisState):
            return self.action(other)
        elif isinstance(other, StateVec):
            return StateVec.from_states(*[
                self.action(b).__rmul__(a) for b, a in other
            ])
        elif isinstance(other, VecOperator):
            return FuncOperator(
                lambda b: self.__mul__(other.__mul__(b)),
                self.is_ketop
            )
        elif isinstance(other, FuncOperator):
            return FuncOperator(
                lambda b: self.__mul__(other.__mul__(b)),
                self.is_ketop
            )
        elif isinstance(other, Number):
            return FuncOperator(
                lambda b: self.__mul__(b).__rmul__(other),
                self.is_ketop
            )

    def __rmul__(self,
            other: BasisState | StateVec | VecOperator | FuncOperator | Number):
        assert isinstance(other, (BasisState, *Number)) \
                or (isinstance(other, StateVec)
                    and not self.is_ketop and not other.is_ket) \
                or (isinstance(other, (VecOperator, FuncOperator))
                    and self.is_ketop == other.is_ketop)
        if isinstance(other, BasisState):
            return self.action(other)
        elif isinstance(other, StateVec):
            return StateVec.from_states(*[
                self.action(b).__rmul__(a) for b, a in other
            ])
        elif isinstance(other, VecOperator):
            return FuncOperator(
                lambda b: self.__rmul__(other.__rmul__(b)),
                self.is_ketop
            )
        elif isinstance(other, FuncOperator):
            return FuncOperator(
                lambda b: self.__rmul__(other.__rmul__(b)),
                self.is_ketop
            )
        elif isinstance(other, Number):
            return FuncOperator(
                lambda b: self.__rmul__(b).__rmul__(other),
                self.is_ketop
            )

    def __truediv__(self, other: Number):
        assert isinstance(other, Number)
        return FuncOperator(
            lambda b: self.action(b) / other,
            self.is_ketop
        )

    def __call__(self, operand: Ket | StateVec):
        assert isinstance(operand, Ket) \
                or (isinstance(operand, StateVec) and operand.is_ket)
        if isinstance(operand, Ket):
            return operand.amp * self.action(operand.state)
        else:
            return StateVec.from_states(*[
                a * self.action(s) for s, a in operand
            ])

    def relabeled(self, label_func):
        return FuncOperator(lambda s: self.action(s).relabeled(label_func))

def _gen_bloch_sphere(labels: list[str, str], is_ket: bool=True) -> pd.Plotter:
    th = np.linspace(0, np.pi, 50)
    ph = np.linspace(0, 2 * np.pi, 50)
    TH, PH = np.meshgrid(th, ph)
    X = np.sin(TH) * np.cos(PH)
    Y = np.sin(TH) * np.sin(PH)
    Z = np.cos(TH)

    my_labels = [f"$\\left| {l} \\right\\rangle$" if is_ket
        else f"$\\left\\langle {l} \\right|$" for l in labels]

    P = pd.Plotter.new_3d(figsize=(2.25, 1.8)) \
        .plot_surface(X, Y, Z, alpha=0.25, color="#f0f0f0") \
        .plot(np.cos(ph), np.sin(ph), np.zeros(ph.shape),
            color="#e8e8e8", linewidth=0.5) \
        .plot(np.sin(ph), np.zeros(ph.shape), np.cos(ph),
            color="#e8e8e8", linewidth=0.5) \
        .plot(np.zeros(ph.shape), np.sin(ph), np.cos(ph),
            color="#e8e8e8", linewidth=0.5) \
        .plot([-1, +1], [0, 0], [0, 0], color="#e8e8e8", linestyle=":",
            linewidth=0.5) \
        .plot([0, 0], [-1, +1], [0, 0], color="#e8e8e8", linestyle=":",
            linewidth=0.5) \
        .plot([0, 0], [0, 0], [-1, +1], color="#e8e8e8", linestyle=":",
            linewidth=0.5) \
        .scatter([0], [0], [0], color="k", s=1) \
        .scatter([1], [0], [0], color="#e0e0e0", s=0.5) \
        .scatter([0], [1], [0], color="#e0e0e0", s=0.5) \
        .scatter([0], [0], [1], color="#e0e0e0", s=0.5) \
        .set_box_aspect((1, 1, 1)) \
        .axis("off") \
        .text(0, 0, +1.1, my_labels[0],
            horizontalalignment="center",
            verticalalignment="bottom",
            fontsize="x-small"
        ) \
        .text(0, 0, -1.1, my_labels[1],
            horizontalalignment="center",
            verticalalignment="top",
            fontsize="x-small"
        )
    return P

def _draw_bloch(alpha: complex, beta: complex, labels: list[str, str],
        P: pd.Plotter) -> pd.Plotter:
    N = np.sqrt(abs(alpha)**2 + abs(beta)**2)
    alpha_ = alpha / N
    beta_ = beta / N
    phi0 = np.arctan2(alpha_.imag, alpha_.real)
    phi1 = np.arctan2(beta_.imag, beta_.real)
    phi = phi1 - phi0
    theta = 2 * np.arccos(abs(alpha_))
    R = min(1, N)

    _draw_bloch_state(R, theta, phi, labels, P)
    return P

def _draw_bloch_state(R: float, theta: float, phi: float,
        labels: list[str, str], P: pd.Plotter):
    x = R * np.sin(theta) * np.cos(phi)
    y = R * np.sin(theta) * np.sin(phi)
    z = R * np.cos(theta)
    P \
        .quiver(
            [0], [0], [0], [x], [y], [z],
            color="k",
            linewidth=0.65,
            arrow_length_ratio=0.1,
            capstyle="round"
        ) \
        .scatter([x], [y], [z], color="C0", s=1)
    a = np.cos(theta / 2)
    b = np.sin(theta / 2)
    X = P.fig.text(0.5, 0.9,
        f"${a:.2f} \\left|{labels[0]}\\right\\rangle \
            {b:+0.2f} e^{{{phi / np.pi:+0.2f} \\pi i}} \\left|{labels[1]}\\right\\rangle$",
        horizontalalignment="center",
        fontsize="xx-small"
    )
    P.outputs.append(X)

def _gen_theta_slider(P: pd.Plotter):
    theta_slider_ax = P.fig.add_axes(
        [0.25, 0.12, 0.375, 0.02],
        facecolor="#e0e0e0"
    )
    theta_slider = wig.Slider(theta_slider_ax,
        label="$\\theta$",
        valmin=0, valmax=1, valinit=0,
        valfmt="$%+.2f \\pi$"
    )
    theta_slider.label.set_fontsize("x-small")
    theta_slider.valtext.set_fontsize("x-small")
    return theta_slider

def _gen_phi_slider(P: pd.Plotter):
    phi_slider_ax = P.fig.add_axes(
        [0.25, 0.05, 0.375, 0.02],
        facecolor="#e0e0e0"
    )
    phi_slider = wig.Slider(phi_slider_ax,
        label="$\\varphi$",
        valmin=-2, valmax=+2, valinit=0,
        valfmt="$%+.2f \\pi$"
    )
    phi_slider.label.set_fontsize("x-small")
    phi_slider.valtext.set_fontsize("x-small")
    return phi_slider

def _gen_reset_button(P: pd.Plotter):
    reset_button_ax = P.fig.add_axes([0.0333, 0.9, 0.15, 0.0667])
    reset_button = wig.Button(reset_button_ax, "Reset",
        color="white", hovercolor="#e0e0e0"
    )
    reset_button.label.set_fontsize("xx-small")
    return reset_button

def draw_bloch_interactive(labels: list[str, str]):
    P = _gen_bloch_sphere(labels, is_ket=True)
    P.fig.subplots_adjust(bottom=0.35)

    theta_slider = _gen_theta_slider(P)
    phi_slider = _gen_phi_slider(P)
    _draw_bloch_state(1.0, 0.0, 0.0, labels, P)
    def sliders_on_changed(val):
        P.outputs.pop().remove()
        P.outputs.pop().remove()
        P.outputs.pop().remove()
        theta = theta_slider.val * np.pi
        phi = phi_slider.val * np.pi
        a = np.cos(theta / 2)
        b = np.sin(theta / 2)
        _draw_bloch_state(1.0, theta, phi, labels, P)
        P.fig.canvas.draw_idle()
    theta_slider.on_changed(sliders_on_changed)
    phi_slider.on_changed(sliders_on_changed)

    reset_button = _gen_reset_button(P)
    def reset_button_on_clicked(mouse_event):
        theta_slider.reset()
        phi_slider.reset()
    reset_button.on_clicked(reset_button_on_clicked)
    P.show()

def _gen_projection_toggle(P: pd.Plotter):
    projection_toggle_ax = P.fig.add_axes([0.0333, 0.8, 0.15, 0.0667])
    projection_toggle = wig.Button(projection_toggle_ax, "Project",
        color="white", hovercolor="#e0e0e0"
    )
    projection_toggle.label.set_fontsize("xx-small")
    return projection_toggle

def _gen_measure_buttons(P: pd.Plotter):
    buttons = list()
    for i in range(3):
        ax = P.fig.add_axes([0.0333 + 0.05 * i, 0.7, 0.05, 0.0667])
        button = wig.Button(ax, "XYZ"[i],
            color="white", hovercolor="#e0e0e0")
        button.label.set_fontsize("xx-small")
        buttons.append(button)
    return buttons

def draw_bloch_interactive_measure(labels: list[str, str]):
    P = _gen_bloch_sphere(labels)
    P.fig.subplots_adjust(bottom=0.35)

    theta_slider = _gen_theta_slider(P)
    phi_slider = _gen_phi_slider(P)
    _draw_bloch_state(1.0, 0.0, 0.0, labels, P)
    def sliders_on_changed(val):
        P.outputs.pop().remove()
        P.outputs.pop().remove()
        P.outputs.pop().remove()
        theta = theta_slider.val * np.pi
        phi = phi_slider.val * np.pi
        a = np.cos(theta / 2)
        b = np.sin(theta / 2)
        _draw_bloch_state(1.0, theta, phi, labels, P)
        P.fig.canvas.draw_idle()
    theta_slider.on_changed(sliders_on_changed)
    phi_slider.on_changed(sliders_on_changed)

    reset_button = _gen_reset_button(P)
    def reset_button_on_clicked(mouse_event):
        theta_slider.reset()
        phi_slider.reset()
        data[:, :] = 0.0
        for i in range(3):
            Pm[i].outputs.pop().remove()
            Pm[i].outputs.pop().remove()
            Pm[i].barh([0, 1], [0.0, 0.0],
                height=0.6, linewidth=0.7, color=f"C{c[i]}", edgecolor="k",
                zorder=3)
            Pm[i].errorbar([0.0, 0.0], [0, 1], xerr=[0.0, 0.0],
                linestyle="", color="k",
                zorder=4)
        Pm.fig.canvas.draw_idle()
    reset_button.on_clicked(reset_button_on_clicked)

    global do_projection
    do_projection = False
    projection_toggle = _gen_projection_toggle(P)
    def projection_button_on_clicked(mouse_event):
        global do_projection
        do_projection = not do_projection
        if do_projection:
            projection_toggle.color="C0"
            projection_toggle.hovercolor="#1c6ba2"
        else:
            projection_toggle.color="white"
            projection_toggle.hovercolor="#e0e0e0"
        P.fig.canvas.draw_idle()
    projection_toggle.on_clicked(projection_button_on_clicked)

    c = [0, 1, 3]
    Pm = pd.Plotter.new(figsize=(1.25, 1.8), nrows=3, sharex=True)
    for i in range(3):
        Pm[i].set_ylabel("XYZ"[i], fontsize="small")
        Pm[i].ggrid()
        Pm[i].grid(False, axis="y", which="minor")
        Pm[i].tick_params(labelsize="x-small")
        Pm[i].set_ylim(-0.5, 1.5)
        Pm[i].set_yticks([0, 1])
        Pm[i].set_yticklabels(["0", "1"])
    Pm[2].set_xlabel("Probability", fontsize="small")
    Pm[2].set_xlim(0, 1.1)
    Pm.tight_layout(pad=0.1, h_pad=0.1, w_pad=0)
    for i in range(3):
        Pm[i].barh([0, 1], [0.0, 0.0],
            height=0.6, linewidth=0.7, color=f"C{c[i]}", edgecolor="k",
            zorder=3)
        Pm[i].errorbar([0.0, 0.0], [0, 1], xerr=[0.0, 0.0],
            linestyle="", color="k", zorder=4)
    data = np.zeros((3, 3))

    measure_x, measure_y, measure_z = _gen_measure_buttons(P)
    def do_measure_x(mouse_event):
        global do_projection
        theta = theta_slider.val * np.pi
        phi = phi_slider.val * np.pi
        p0 = (1 + np.sin(theta) * np.cos(phi)) / 2
        B = 1 - int(random.random() < p0)
        data[0][B] = (data[0][2] * data[0][B] + 1) / (data[0][2] + 1)
        data[0][1 - B] = data[0][2] * data[0][1 - B] / (data[0][2] + 1)
        data[0][2] += 1
        Pm[0].outputs.pop().remove()
        Pm[0].outputs.pop().remove()
        Pm[0].barh([0, 1], data[0][:2],
            height=0.6, linewidth=0.7, color="C0", edgecolor="k", zorder=3)
        Pm[0].errorbar(data[0][:2], [0, 1],
            xerr=np.sqrt(data[0][:2]) / data[0][2],
            color="k", linestyle="", zorder=4)
        if do_projection:
            theta_slider.set_val(0.5)
            phi_slider.set_val(1 if B else 0)
        Pm.fig.canvas.draw_idle()
    measure_x.on_clicked(do_measure_x)

    def do_measure_y(mouse_event):
        global do_projection
        theta = theta_slider.val * np.pi
        phi = phi_slider.val * np.pi
        p0 = (1 + np.sin(theta) * np.sin(phi)) / 2
        r = random.random()
        print(r, p0)
        B = 1 - int(r < p0)
        print(B)
        data[1][B] = (data[1][2] * data[1][B] + 1) / (data[1][2] + 1)
        data[1][1 - B] = data[1][2] * data[1][1 - B] / (data[1][2] + 1)
        data[1][2] += 1
        Pm[1].outputs.pop().remove()
        Pm[1].outputs.pop().remove()
        Pm[1].barh([0, 1], data[1][:2],
            height=0.6, linewidth=0.7, color="C1", edgecolor="k", zorder=3)
        Pm[1].errorbar(data[1][:2], [0, 1],
            xerr=np.sqrt(data[1][:2]) / data[1][2],
            color="k", linestyle="", zorder=4)
        if do_projection:
            theta_slider.set_val(0.5)
            phi_slider.set_val(-0.5 if B else 0.5)
        Pm.fig.canvas.draw_idle()
    measure_y.on_clicked(do_measure_y)

    def do_measure_z(mouse_event):
        global do_projection
        theta = theta_slider.val * np.pi
        phi = phi_slider.val * np.pi
        p0 = np.cos(theta / 2)**2
        B = 1 - int(random.random() < p0)
        data[2][B] = (data[2][2] * data[2][B] + 1) / (data[2][2] + 1)
        data[2][1 - B] = data[2][2] * data[2][1 - B] / (data[2][2] + 1)
        data[2][2] += 1
        Pm[2].outputs.pop().remove()
        Pm[2].outputs.pop().remove()
        Pm[2].barh([0, 1], data[2][:2],
            height=0.6, linewidth=0.7, color="C3", edgecolor="k", zorder=3)
        Pm[2].errorbar(data[2][:2], [0, 1],
            xerr=np.sqrt(data[2][:2]) / data[2][2],
            color="k", linestyle="", zorder=4)
        if do_projection:
            theta_slider.set_val(B)
        Pm.fig.canvas.draw_idle()
    measure_z.on_clicked(do_measure_z)

    pp.show()

def draw_state_interactive(labels: list[str, str]):
    draw_bloch_interactive(labels)

