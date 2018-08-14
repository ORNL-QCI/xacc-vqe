import pyxaccvqe as vqe
from pyxaccvqe import qpu, PauliOperator
import pyxacc as xacc
import numpy as np

xacc.Initialize()

ham = PauliOperator(5.906709445) + \
    PauliOperator({0:'X',1:'X'}, -2.1433) + \
    PauliOperator({0:'Y',1:'Y'}, -2.1433) + \
    PauliOperator({0:'Z'}, .21829) + \
    PauliOperator({1:'Z'}, -6.125)

@qpu.vqe(accelerator='tnqvm', observable=ham, optimizer='scipy-nelder-mead', opt_params={'tol':1e-2})
def foo(initialTheta):
   X(0)
   Ry(initialTheta, 1)
   CNOT(1, 0)
   return
   
result = foo()
print(result.energy, result.angles)

# Get the expectation values for the Z0 term
@qpu(accelerator='tnqvm')
def foo(theta):
   X(0)
   Ry(theta, 1)
   CNOT(1, 0)
   Measure(0,0)
   return

expVals = [foo(t).getExpectationValueZ() for t in np.linspace(-np.pi,np.pi,10)]
print (expVals)

xacc.Finalize()