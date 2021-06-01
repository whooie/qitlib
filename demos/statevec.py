import numpy as np
import lib

# define a basis - everything in qitlib assumes we're only working in one basis
B = lib.Basis.from_labels("+z", "-z")

# create some StateVecs
# creation methods are flexible: can use the StateVec constructor method, a math
# expression, or directly from the BasisStates
z_p = lib.StateVec({B[0]: 1}) # |+z>
z_m = B[1].to_statevec() # |-z>
y_p = lib.StateVec({B[0]: 1, B[1]: 1j}).normalized() # |+y>
y_m = (z_p - 1j * z_m) / np.sqrt(2) # |-y>

# create some operators
S_z = lib.VecOperator({B[0]: z_p / 2, B[1]: -z_m / 2})
S_x = lib.VecOperator({B[0]: z_m / 2, B[1]: z_p / 2})
S_y = lib.VecOperator({B[0]: 1j * z_m / 2, B[1]: -1j * z_p / 2})

# StateVecs print out in Dirac notation
print("PRINTING")
print("|+z> =", z_p)
print("|+y> =", y_p)
print("<-y| =", y_m.hc())

print("")

# do some math
print("MATH")
print("(|+y> + |-y>) / sqrt(2) =", (y_p + y_m) / np.sqrt(2))
print("i|+y> =", 1j * y_p)
print("|+z> - |+z> =", z_p - z_p)
print("<+z|+y> =", z_p.hc() * y_p)
print("<-z|-y> =", z_m.hc() * y_m)

print("")

# do some math
print("OPERATOR MATH")
print("S_z =", S_z)
print("S_z * S_y =", S_z * S_y)
print("S_z |+z> =", S_z * z_p)
print("S_z |+y> =", S_z * y_p)
# Python's order of operations makes expectation values a bit awkward
print("<-z| S_z |-z> =", z_m.hc() * (S_z * z_m))
print("<+y| S_y |+y> =", y_p.hc() * (S_y * y_p))
print("<-z| S_x |+z> =", z_m.hc() * S_x.hc() * z_p)

