#pragma vqe-coefficient 15.5317
__qpu__ idTerm(AcceleratorBuffer b) {
}
#pragma vqe-coefficient -2.1433
__qpu__ x0x1(AcceleratorBuffer b) {
H 0
H 1
MEASURE 0 [0]
MEASURE 1 [1]
}
#pragma vqe-coefficient -2.1433
__qpu__ y0y1(AcceleratorBuffer b) {
RX(1.57078) 0
RX(1.57078) 1
MEASURE 0 [0]
MEASURE 1 [1]
}
#pragma vqe-coefficient .21829
__qpu__ z0(AcceleratorBuffer b) {
MEASURE 0 [0]
}
#pragma vqe-coefficient -6.125
__qpu__ z1(AcceleratorBuffer b) {
MEASURE 1 [0]
}
#pragma vqe-coefficient -9.625
__qpu__ z2(AcceleratorBuffer b) {
MEASURE 2 [0]
}
#pragma vqe-coefficient -3.91312
__qpu__ x1x2(AcceleratorBuffer b) {
H 1
H 2
MEASURE 1 [0]
MEASURE 2 [1]
}
#pragma vqe-coefficient -3.91312
__qpu__ y1y2(AcceleratorBuffer b) {
RX(1.57078) 1
RX(1.57078) 2
MEASURE 1 [0]
MEASURE 2 [1]
}
