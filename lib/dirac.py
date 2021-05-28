from __future__ import annotations
import numpy as np
import scipy.sparse as sp
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
        self.components = components
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
                    _components[b] = _components[b] + a \
                            if b in _components.keys() else a
                    if _components[b] == 0:
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
            return " + ".join(f"({a:g})" \
                    + ("|" if self.is_ket else "<") \
                    + f"{b.label}" \
                    + (">" if self.is_ket else "|")
                for b, a in self if a != 0
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

class VecOperator:
    """
    Quantum operator represented as a collection of `BasisState`-to-`StateVec`
    mappings. The action of the operator is internally stored as
    dict[BasisState: StateVec] (and is indexed as such). Interacts only with
    other `VecOperator`s, `FuncOperator`s, and `StateVec`s.
    """
    def __init__(self, action: dict[BasisState, StateVec], is_ketop: bool=True,
            default_zero: bool=True, scalar: complex=1.0 + 0.0j):
        self.action = action
        self.is_ketop = is_ketop
        self.default_zero = default_zero
        self.X = complex(scalar)

    @staticmethod
    def identity(is_ketop=True):
        return VecOperator(
            dict(),
            is_ketop,
            default_zero=False,
            scalar=1.0 + 0.0j
        )

    def adjoint(self):
        return VecOperator(
            {b: self.action[s].conjugate() for b in self.keys()},
            not self.is_ketop,
            self.default_zero,
            self.X
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
        return self.X * self.action.get(key, default)

    def __getitem__(self, b: BasisState):
        assert isinstance(b, BasisState)
        return self.X * self.action.get(
            b,
            StateVec.zero(self.is_ketop) if self.default_zero \
                    else StateVec({b: 1.0 + 0.0j}, self.is_ketop)
        )

    def __neg__(self):
        return VecOperator(
            {b: s for b, s in self},
            self.is_ketop,
            self.default_zero,
            -self.X
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
                {b: self[b] + other.X * other[b] for b in all_keys},
                self.is_ketop,
                default_zero=True,
                scalar=1.0 + 0.0j
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
                self.default_zero * other.default_zero,
                self.X * other.X
            )
        elif isinstance(other, Number):
            return VecOperator(
                self.action,
                self.is_ketop,
                self.default_zero,
                self.X * other
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
                other.default_zero * self.default_zero,
                other.X * self.X
            )
        elif isinstance(other, Number):
            return VecOperator(
                {b: self[b].__rmul__(other) for b in self.keys()},
                self.is_ketop,
                self.default_zero,
                other * self.X
            )
        else:
            return other.__mul__(self)

    def __truediv__(self, other: Number):
        assert isinstance(other, Number)
        return VecOperator(
            self.action,
            self.is_ketop,
            self.default_zero,
            self.X / other
        )

    def __iter__(self):
        return iter(self.action.items())

    def relabeled(self, label_func):
        return VecOperator(
            {b.relabeled(label_func): s.relabeled(label_func) for b, s in self},
            self.is_ketop,
            self.default_zero,
            self.X
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
            scalar=1.0 + 0.0j
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
        action(BasisState) -> StateVec(is_ket=True)
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
            (lambda b: StateVec({b: 1.0 + 0.0j})) if is_ketop \
                    else (lambda b: StateVec({b: 1.0 + 0.0j}, is_ket=False)),
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
            scalar=1.0 + 0.0j
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

