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
the standard arithmetic operations. Once a `StateVec` or `StateArr` has been
created, it can easily be converted between ket and bra forms with the standard
`.conjugate()` or shortened `.hc()` methods.

In this library, however, ket-bras (e.g. `|psi><phi|`) are not allowed, and will
raise an exception. Additionally, both types provide methods to calculate their
normalization with `.norm()` (and hence compute an equivalent, properly
normalized version of the state with `.normalized()`) as well as sample their
corresponding probability distributions with `.measure()` and
`.measure_single()`.

Choosing which representation to use comes down mostly to details regarding
computation speed and memory usage that are a bit outside the scope of this
program. Both types can be easily convert to and from each other with
`.to_statevec()` and `to_statearr()`, so feel free to have fun with them!

## 3 Operators
Like `StateVec`s and `StateArr`s, quantum-mechanical operators are also
represented in different ways. Here, there are three: `VecOperator`,
`MatOperator`, and `FuncOperator`.

### 3.1 Different implementations of operators
The difference between `VecOperator` and `MatOperator` is almost exactly the
same as that between `StateVec` and `StateArr`: the first is represented
internally by a `dict` mapping `BasisState`s to `StateVec`s, while the second
uses a matrix. `FuncOperator` is a bit different, though, so we'll come back to
it in a bit.

Both kinds of operators obey the regular mathematical addition, subtraction,
multiplication, and division rules with numbers and other operators of the same
kind. But as their names suggest, `VecOperator`s can only multiply `StateVec`s
and `MatOperator`s can only multiply `StateArr`s, although like states, they can
be easily converted between forms with `.to_matoperator()` or
`.to_vecoperator()`. Both kinds of operators can also be converted to their
adjoint forms with the `.adjoint()` or `.hc()` methods.

The action performed by these operators is only defined in a single direction.
That is, adjoint operators (which have field `.is_ketop = False`) can only
multiply bras (states with field `.is_ket = False`) from the right, and
non-adjoint operators (`.is_ketop = True`) can only multiply kets (`.is_ket =
True`) from the left. One gotcha to keep in mind is that Python always prefers
to apply operations from the left regardless of the types involved, so liberal
use of parentheses is encouraged when multiplying operators and states.

As with the state object types, the choice of which form of operator to use
mostly comes down to preference. One limitation of `MatOperator`s and
`VecOperator`s, however, is that their action may only be defined over a finite
basis. `FuncOperator`s provide a solution to this by having their actions be
represented internally as a Python function, rather than a `dict` or an array.
`FuncOperator`s can only multiply `StateVec`s, and their functions must have the
type signature
```Python
action(BasisState) -> StateVec(is_ket=self.is_ketop)
```
(that is, the function must accept only one argument of type `BasisState`, and
return only one `StateVec` whose `.is_ket` field must have the same value as the
`.is_ketop` field of the operator).

## 4 Circuits


### 4.1 Atomic and non-atomic gates

### 4.2 Reading and writing circuit files

### 4.3 Circuit output

### 4.4 Drawing circuits

[1]: https://qiskit.org/

