from lib.dirac import *
from lib.qcomp import *
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pp
import random
from fractions import Fraction

def is_prime(n):
    """
    Return True if `n` is a prime number and False otherwise.
    """
    return not (n < 2 or any(n % x == 0 for x in range(2, int(np.sqrt(n)) + 1)))

def is_pow(n):
    """
    Return True if `n` is a perfect power of a number and False otherwise.
    """
    return any(pow(n, 1/a) % 1 == 0 for a in range(2, int(np.log2(n)) + 1))

def gcd(a, b):
    """
    Return the greatest common divisor of `a` and `b` using Euclid's algorithm.
    """
    q = b // a
    r = b % a
    if r == 0:
        return a
    else:
        return gcd(r, a)

def order(x, N, maxiters=10000):
    """
    Return the number-theoretical order of `x`, mod `N`. (That is, the smallest
    power `r` such that
        x^r == 1 mod N
    """
    for r in range(1, maxiters + 1):
        if pow(x, r, N) == 1:
            return r

def shor_classical(N, maxiters=10000, xr: bool=False, allow_lucky: bool=False,
        print_flag: bool=False):
    """
    Implement Shor's algorithm using only classical operations.
    """
    assert not is_prime(N)
    assert not is_pow(N)
    assert N % 2 != 0
    for k in range(maxiters):
        if print_flag:
            print(f"\riter = {k + 1} ", end="", flush=True)
        x = int((np.sqrt(N) - 1) * random.random() + 2)
        if print_flag:
            print(f" x = {x} ", end="", flush=True)
        if (a := gcd(x, N)) != 1:
            if allow_lucky:
                if xr:
                    return a, N // a, -1, -1
                else:
                    return a, N // a
            else:
                continue
        r = order(x, N)
        if print_flag:
            print(f" r = {r}          ", end="", flush=True)
        if r % 2 == 1:
            continue
        a = (x**(r // 2) - 1) % N
        b = (x**(r // 2) + 1) % N
        if 0 in {a, b} or 1 in {a, b}:
            continue
        A = gcd(a, N)
        B = gcd(b, N)
        if print_flag:
            print("")
        if xr:
            return A, B, x, r
        else:
            return A, B

def eigenvalue_period(U, N, maxiters=10000):
    e = (np.log(la.eigvals(U)) / (2 * np.pi * 1j)).real
    r = 1
    for j in range(maxiters):
        r = Fraction(random.choice(e)).limit_denominator(N).denominator
        if r != 1:
            break
    return r

def shor_semiclassical(N, maxiters=10000, xr: bool=False,
        allow_lucky: bool=False, print_flag: bool=False):
    """
    Implement Shor's algorithm semiclassically: use all classical operations,
    except for the determination of the order r, which is computed using the
    matrix representation of the quantum-mechanical xymodN operator.
    """
    assert not is_prime(N)
    assert not is_pow(N)
    assert N % 2 != 0
    for k in range(maxiters):
        if print_flag:
            print(f"\riter = {k + 1} ", end="", flush=True)
        x = int((np.sqrt(N) - 1) * random.random() + 2)
        if print_flag:
            print(f" x = {x} ", end="", flush=True)
        if (a := gcd(x, N)) != 1:
            if allow_lucky:
                if xr:
                    return a, N // a, -1, -1
                else:
                    return a, N // a
            else:
                continue
        n = int(np.ceil(np.log2(N)))
        U = xymodN_mat(n, 0, n, x, N)
        e = (np.log(la.eigvals(U)) / (2 * np.pi * 1j)).real
        r = eigenvalue_period(U, N, maxiters)
        if print_flag:
            print(f" r = {r}          ", end="", flush=True)
        if r % 2 == 1:
            continue
        a = (x**(r // 2) - 1) % N
        b = (x**(r // 2) + 1) % N
        if 0 in {a, b} or 1 in {a, b}:
            continue
        A = gcd(a, N)
        B = gcd(b, N)
        if A * B != N:
            continue
        if print_flag:
            print("")
        if xr:
            return A, B, x, r
        else:
            return A, B

def shor_quantum(N, maxiters=10000, xr: bool=False, allow_lucky: bool=False,
        print_flag: bool=False):
    """
    Implement Shor's algorithm using a simulated quantum circuit to determine
    the order r.
    """
    assert not is_prime(N)
    assert not is_pow(N)
    assert N % 2 != 0
    for k in range(maxiters):
        if print_flag:
            print(f"\riter = {k + 1} ", end="", flush=True)
        x = int((np.sqrt(N) - 1) * random.random() + 2)
        if print_flag:
            print(f" x = {x} ", end="", flush=True)
        if (a := gcd(x, N)) != 1:
            if allow_lucky:
                if xr:
                    return a, N // a, -1, -1
                else:
                    return a, N // a
            else:
                continue
        n = int(np.ceil(np.log2(N)))
        nwires = 3*n + 1
        basis = qubit_basis(nwires)
        qft_wires = nwires - n
        U = xymodN_mat(n, 0, n, x, N)
        initial_state = StateArr(
            np.append(
                la.eig(U)[1][:, random.randrange(0, 2**n)],
                np.zeros(2**(nwires) - 2**n),
            ),
            basis,
            is_ket=True
        )
        C = Circuit(
            nwires=nwires,
            elements=[("H", i) for i in range(qft_wires)] \
                    + [("CFUNC", qft_wires - 1 - j, qft_wires, n, "xyModN", x**(2**j), N)
                        for j in range(qft_wires)] \
                    + iqft_elements(0, qft_wires - 1) \
                    + reverse_elements(0, qft_wires - 1),
            initial_state=initial_state,
            measure=False
        )
        output = C.run(kind="func", stepwise=True)
        state = output.get_statevec().most_probable()
        if print_flag:
            print(f" state = {state} ", end="", flush=True)
        theta = b2i(state.label[:qft_wires]) / 2**qft_wires
        r = Fraction(theta).limit_denominator(N).denominator
        if print_flag:
            print(f" r = {r}   ", end="", flush=True)
        if r % 2 == 1:
            continue
        a = (x**(r // 2) - 1) % N
        b = (x**(r // 2) + 1) % N
        if 0 in {a, b} or 1 in {a, b}:
            continue
        A = gcd(a, N)
        B = gcd(b, N)
        if A * B != N:
            continue
        if print_flag:
            print("")
        if xr:
            return A, B, x, r
        else:
            return A, B

