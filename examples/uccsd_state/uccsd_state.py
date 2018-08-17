import xacc
import numpy as np

xacc.Initialize() 
tnqvm = xacc.getAccelerator('tnqvm')

nQubits = 4
nElectrons = 2

uccsdGen = xacc.getIRGenerator("uccsd")
f = uccsdGen.generate({'n-electrons':nElectrons, 'n-qubits':nQubits})

f = f.eval([3.14,-1.51])
state = xacc.gate.getState(tnqvm, f)

print(state)