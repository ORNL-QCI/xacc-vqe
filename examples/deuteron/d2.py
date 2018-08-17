import sys
import xacc
import xaccvqe
from xaccvqe import PauliOperator

xacc.Initialize(['--compiler','quil'])

ham = PauliOperator(5.906709445) + \
        PauliOperator({0:'X',1:'X'}, -2.1433) + \
        PauliOperator({0:'Y',1:'Y'}, -2.1433) + \
        PauliOperator({0:'Z'}, .21829) + \
        PauliOperator({1:'Z'}, -6.125)

print('\nH_{2x2} = ', ham)

@xaccvqe.qpu.vqe(accelerator='tnqvm', observable=ham)
def ansatz(t0):
    X(0)
    Ry(t0, 1)
    CNOT(1, 0)


vqeResult = ansatz() #vqe.execute(ham, **{'task':'vqe', 'ansatz':ansatz})
print('(Optimal Angle, Energy) = (', vqeResult.angles, ',', vqeResult.energy, ')')
print('Number of QPU Calls = ', vqeResult.nQpuCalls)
print('Number of VQE Iterations = ', vqeResult.vqeIterations)

xacc.Finalize()
