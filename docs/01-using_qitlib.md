# Using `qitlib`

`qitlib` is a small library I wrote for simulating simple quantum-computational
circuits. This document intends to describe the high-level features of the
library as they pertain to pedagogical usage. This library is far from
full-featured; a richer (if somewhat less accessible) experience can be obtained
with e.g. [Qiskit][1]. Basic knowledge of Python is assumed.

More specific information regarding things like exact usage and function strings
can be found via the docstrings in the appropriate library files. A few examples
can be found in the `demos` folder.

## 1 Bases
To do any quantum computation, you first need a basis. From a computational
perspective, a basis defines the space of all the results one could possibly
arrive at using any given sequence of operations. More simply, a basis is the
set of objects in terms of which everything else in our computations can be
written.

Here in `qitlib`, a `Basis` is a list-like object whose contents are the
`BasisState`s that define the computational space. In particular, `qitlib`
assumes that all of our operations are done using only a single basis of
orthonormal states. That is, all computations will be done in terms of only a
single set of basis states, which all multiply each other to 0 and themselves to
1.

Essentially, `Basis` objects act like `list`s that only hold `BasisState`s,
while the `BasisState`s hold the labels associated with each state.

## 2 Quantum states
`qitlib` provides to ways to represent quantum states. The first is as a sum of
`BasisState`s, while the second is as a list of coordinates corresponding to the
elements of a given `Basis`.

### 2.1 `StateVec`s vs `StateArr`s
The internal representation of the quantum state is where these two classes
differ. Under the hood, `StateVec`s store their data in a `dict` that maps
`BasisState`s to their corresponding complex amplitudes. On the other hand,
`StateArr`s store their data as an array of coordinates along with a `Basis`
that tells them which coordinate goes with which `BasisState`.

Both object types act like their underlying data structures (i.e. `StateVec`s
can be fed `BasisState` keys and `StateArr`s can be indexed like arrays, and
both can be iterated over in the same way using `for x in y` syntax), along with
the standard arithmetic operations. In this case, however, ket-bras
(`|psi><phi|`) are not allowed, and will raise an exception. Additionally, both
types provide methods to calculate their normalization (and hence compute an
equivalent, properly normalized version of the state) as well as sample their
corresponding probability distributions.

Choosing which representation to use comes down to mostly irrelevant details
about computation speed and memory usage.

## 3 Operators
Like `StateVec`s and `StateArr`s, quantum-mechanical operators are also
represented in different ways.

### 3.1 Different implementations of operators

## 4 Circuits

### 4.1 Atomic and non-atomic gates

### 4.2 Reading circuits

### 4.3 Circuit output

### 4.4 Drawing circuits

[1]: https://qiskit.org/

