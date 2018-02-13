import numpy as np
import pyxacc as xacc
from pyxacc import InstructionParameter
import pyxaccvqe as vqe
from pyxaccvqe import PauliOperator
from scipy.optimize import minimize

xacc.Initialize()

# Construct the First Quantized 2x2 and 3x3 Hamiltonians
hamiltonian3x3 = PauliOperator(7.7658547225) + PauliOperator({0:'X'}, -2.143303525) + \
                PauliOperator({0:'X', 1:'X'}, -3.91311896) + PauliOperator({0:'X', 1:'Z'}, -2.143303525) + \
                PauliOperator({0:'Y',1:'Y'}, -3.91311896) + PauliOperator({0:'Z'}, 1.6408547224999999) + \
                PauliOperator({0:'Z',1:'Z'}, -7.9841452775) + PauliOperator({1:'Z'}, -1.8591452775000001)
            
print('\nH_{3x3} = ', hamiltonian3x3)

# Create the 3x3 Ansatz
ansatz3x3 = xacc.gate.GateFunction('statePrep', ['theta0','theta1'])
ansatz3x3.add(xacc.gate.create('Ry',[0],['theta0']))
ansatz3x3.add(xacc.gate.create('Ry',[1],['theta1']))
print('3x3 Ansatz = \n', ansatz3x3.toString('q'))

qpu = xacc.getAccelerator('tnqvm')

def energy(params):
    ir = hamiltonian3x3.toXACCIR()
    kernels = ir.getKernels()
    qubits = qpu.createBuffer('q',2)
    energy = 0.0
    for k in kernels:
        val = 0
        coeff = k.getParameter(0)
        if k.nInstructions() > 0:
            evaledAnsatz = vqe.AnsatzEvaluator.evaluate(ansatz3x3, 2, np.asarray(params))
            k.insertInstruction(0, evaledAnsatz)
            qpu.execute(qubits, k)
            exp = qubits.getExpectationValueZ()
            qubits.resetBuffer()
            val = coeff * exp
        else:
            val = coeff
        energy += val
    return energy.real

print('XACC Diagonalize: ', vqe.execute(hamiltonian3x3, **{}).energy) #'task':'vqe', 'ansatz':ansatz3x3}).energy)
print('XACC Nelder-Mead: ', vqe.execute(hamiltonian3x3, **{'task':'vqe', 'vqe-params':'.7,.2', 'ansatz':ansatz3x3}).energy)
print('XACC SciPy Minimze:\n', minimize(energy, [0,0]))



from bayes_opt import BayesianOptimization

def neg_energy(x, y):
    return -1*energy([x,y])

bo = BayesianOptimization(neg_energy,
                          {'x': (0, 1), 'y': (0, 1)})

bo.maximize(init_points=10, n_iter=10, kappa=2)

#{'max': {'max_val': 2.035286440109728, 'max_params': {'x': 0.6745710475274774, 'y': 0.2651970328890534}},
#print(bo.res)

print('\nMin Energy = ', -1.0*bo.res['max']['max_val'], ' at ', bo.res['max']['max_params'])
