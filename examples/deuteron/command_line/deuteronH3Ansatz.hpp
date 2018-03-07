__qpu__ ansatz(AcceleratorBuffer b, double t0, double t1) {
X 0
RY(t0) 2
CNOT 2 0
RY(t1) 1
CNOT 0 1
RY(-t1) 1
CNOT 0 1
CNOT 1 0
}

