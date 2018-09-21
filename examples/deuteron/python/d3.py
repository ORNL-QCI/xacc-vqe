import xacc
import xaccvqe
from xaccvqe import PauliOperator

xacc.Initialize(['--compiler','quil'])

ibm = xacc.getAccelerator('local-ibm')
tnqvm = xacc.getAccelerator('tnqvm')
rigetti = xacc.getAccelerator('rigetti')

#buffer = ibm.createBuffer('q', 2)
buffer = tnqvm.createBuffer('q', 2)

h2 = PauliOperator(5.906709445) + \
        PauliOperator({0:'X',1:'X'}, -2.1433) + \
        PauliOperator({0:'Y',1:'Y'}, -2.1433) + \
        PauliOperator({0:'Z'}, .21829) + \
        PauliOperator({1:'Z'}, -6.125)

h3 = h2 + PauliOperator(9.625) + \
    PauliOperator({1:'X',2:'X'}, -3.913119) + \
    PauliOperator({1:'Y',2:'Y'}, -3.913119) + \
    PauliOperator({2:'Z'}, -9.625)


@xaccvqe.qpu.vqe(accelerator=tnqvm, observable=h3)
def ansatz(buffer, t0, t1):
    X(0)
    Ry(t0, 2)
    CNOT(2, 0)
    Ry(t1, 1)
    CNOT(0,1)
    Ry(-t1, 1)
    CNOT(0,1)
    CNOT(1,0)

# Run VQE with given ansatz kernel
initAngle = .2

ansatz(buffer, initAngle, initAngle)

print(buffer.listExtraInfoKeys())
print('Energy = ', buffer.getInformation('vqe-energy'))
print('Opt Angles = ', buffer.getInformation('vqe-angles'))
