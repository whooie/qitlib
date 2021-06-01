import numpy as np
import matplotlib.pyplot as pp
import lib

nwires = 2
B = lib.qubit_basis(nwires) # [ |00>, |01>, |10>, |11> ]

C = lib.Circuit(
    nwires=nwires,
    elements=[
        ("H", 0), # hadamard the zeroth qubit
        ("P", 0, np.pi), # conditional pi-phase shift on the zeroth qubit
        ("NOT", 1), # NOT the first qubit
        ("H", 0), # hadamard the zeroth qubit
    ],
    initial_state=B[0].to_statearr(B), # initial state is |00>
    measure=True # generate 2^nwires * 1000 measurements of the final state
)

# write the circuit to a file; this file can be read by another program later
C.to_file("output/test_circuit.txt")

# can "compile" the circuit to a single unitary
U = C.generate(kind="mat")
print(U.action)

# run the circuit gate-by-gate with VecOperators
output = C.run(stepwise=True, kind="vec")

# stick the measurements in a histogram
output.get_hist(density=True, fractional=False) \
    .savefig("output/test_circuit_output.png") \
    .close()

