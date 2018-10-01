import sys
import xacc
import xaccvqe
from xaccvqe import PauliOperator
import numpy as np

xacc.Initialize()

tnqvm = xacc.getAccelerator('tnqvm')
buffer = tnqvm.createBuffer('q', 2)

ham = PauliOperator(5.906709445) + \
        PauliOperator({0:'X',1:'X'}, -2.1433) + \
        PauliOperator({0:'Y',1:'Y'}, -2.1433) + \
        PauliOperator({0:'Z'}, .21829) + \
        PauliOperator({1:'Z'}, -6.125)

xacc.setOption('itensor-svd-cutoff','1e-16')

# Hardware Efficient Ansatz xacc kernel
@xaccvqe.qpu.energy(accelerator=tnqvm, observable=ham)
def ansatz(buffer, *args):
    xacc(hwe,layers=2,n_qubits=2,connectivity='[[0,1]]')

# XACC Kernels can display number of required nParameters
# and be persisted to qasm string
print(ansatz.nParameters()) 
print(ansatz.getFunction().toString('q'))

# Generate an initial random set of vqe params 
# of the correct number of parameters and execute
ansatz(buffer, *np.random.uniform(low=-np.pi,high=np.pi, size=(ansatz.nParameters(),)))

print(buffer.listExtraInfoKeys())
print('Energy = ', buffer.getInformation('vqe-energy'))
print('Opt Angles = ', buffer.getInformation('vqe-angles'))

xacc.Finalize()
