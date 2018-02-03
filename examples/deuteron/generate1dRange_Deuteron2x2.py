import sys
sys.path.append('/usr/local/xacc/lib/python')
import pyxacc as xacc
import pyxaccvqe as vqe 
import numpy as np
from pyxaccvqe import PauliOperator
from pyxacc import InstructionParameter

# Construct the First Quantized 2x2 and 3x3 Hamiltonians
hamiltonian2x2Id = PauliOperator(5.906709445) 
hamiltonian2x2 = PauliOperator({0:'X'}, -4.28660705) + \
                 PauliOperator({0:'Z'}, -6.343290555)

print('H_{2x2} = ', hamiltonian2x2)

# Create the 2x2 Ansatz
ansatz2x2 = xacc.gate.GateFunction('statePrep', [InstructionParameter('theta')])
ansatz2x2.add(xacc.gate.create('Ry',[0],[InstructionParameter('theta')]))
print('\n2x2 Ansatz = \n', ansatz2x2.toString('q'))

range_ = np.linspace(-np.pi,np.pi,50)
print(range_)

vqeOpts = {'error-mitigation':['correct-readout-errors'], 'n-qubits':2, 'range':[range_], 'ansatz':ansatz2x2}

kernels, irVec = vqe.generateKernels(hamiltonian2x2, **vqeOpts)

print('NKernels = ', len(kernels))

qpu = xacc.getAccelerator('tnqvm')
qubits = qpu.createBuffer('qubits',2)

results = kernels.execute(qubits, [])
print('N Results = ', len(results))

es = [r.getExpectationValueZ() for r in results]
print(es)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

energies = [c[0]*-6.3432290555 + c[1]*-4.28660705 + 5.906709445 for c in chunks(es,2)]

print(energies)

