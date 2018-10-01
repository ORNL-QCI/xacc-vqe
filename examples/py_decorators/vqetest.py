import xaccvqe as vqe
from xaccvqe import qpu, PauliOperator
import xacc
import numpy as np

xacc.Initialize()
tnqvm = xacc.getAccelerator('tnqvm')
buffer = tnqvm.createBuffer('q',2)

ham = PauliOperator(5.906709445) + \
    PauliOperator({0:'X',1:'X'}, -2.1433) + \
    PauliOperator({0:'Y',1:'Y'}, -2.1433) + \
    PauliOperator({0:'Z'}, .21829) + \
    PauliOperator({1:'Z'}, -6.125)

@qpu.vqe(accelerator=tnqvm, observable=ham, optimizer='vqe-bayesopt', opt_params={'tol':1e-2})
def ansatz(buffer, initialTheta):
   X(0)
   Ry(initialTheta, 1)
   CNOT(1, 0)
   
ansatz(buffer, .5)

print(buffer.getInformation('vqe-energy'), buffer.getInformation('vqe-angles'))

xacc.Finalize()