import numpy as np
import lib

# define a basis - everything in qitlib assumes we're only working in one basis
B = lib.Basis.from_labels("+z", "-z")

# create some StateArrs
# creation methods are less convenient: best way is to convert a StateVec
# or BasisState, or use a math expression
z_p = lib.StateArr(np.array([1, 0]), B) # |+z>
z_m = B[1].to_statearr(B) # |-z>
y_p = (z_p + 1j * z_m) / np.sqrt(2) # |+y>
y_m = (z_p - 1j * z_m).normalized() # |-y>

# create some regular spin operators
S_z = lib.MatOperator(np.array([[1, 0], [0, -1]]) / 2, B)
S_x = lib.MatOperator(np.array([[0, 1], [1, 0]]) / 2, B)
S_y = lib.MatOperator(np.array([[0, -1j], [1j, 0]]) / 2, B)

# StateVecs print out in Dirac notation
print("PRINTING")
print("|+z> =")
print(z_p)
print("|+y> =")
print(y_p)
print("<-y| =")
print(y_m.hc())

print("")

# do some math
print("MATH")
print("(|+y> + |-y>) / sqrt(2) =")
print((y_p + y_m) / np.sqrt(2))
print("i|+y> =")
print(1j * y_p)
print("|+z> - |+z> =")
print(z_p - z_p)
print("<+z|+y> =", z_p.hc() * y_p)
print("<-z|-y> =", z_m.hc() * y_m)

print("")

# do some math
print("OPERATOR MATH")
print("S_z =")
print(S_z)
print("S_z * S_y =")
print(S_z * S_y)
print("S_z |+z> =")
print(S_z * z_p)
print("S_x |+z> =")
print(S_x * z_p)
print("S_z |+y> =")
print(S_z * y_p)
# Python's order of operations makes expectation values a bit awkward
print("<-z| S_z |-z> =", z_m.hc() * (S_z * z_m))
print("<+y| S_y |+y> =", y_p.hc() * (S_y * y_p))
print("<-z| S_x |+z> =", z_m.hc() * S_x.hc() * z_p)

