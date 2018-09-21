import sys
import xacc
import xaccvqe
from xaccvqe import PauliOperator

xacc.Initialize(['--compiler','quil'])

ibm = xacc.getAccelerator('local-ibm')
tnqvm = xacc.getAccelerator('tnqvm')

#buffer = ibm.createBuffer('q', 2)
buffer = tnqvm.createBuffer('q', 2)

ham = PauliOperator(5.906709445) + \
        PauliOperator({0:'X',1:'X'}, -2.1433) + \
        PauliOperator({0:'Y',1:'Y'}, -2.1433) + \
        PauliOperator({0:'Z'}, .21829) + \
        PauliOperator({1:'Z'}, -6.125)

print('\nH_{2x2} = ', ham)

xacc.setOption('ibm-shots','8192')
#xacc.setOption('vqe-backend','vqe-bayesopt')
#xacc.setOption('bo-n-iter','20')

def ansatz2(buffer, t0):
    X(0)
    Ry(t0, 1)
    CNOT(1, 0)

@xaccvqe.qpu.vqe(accelerator=tnqvm, observable=ham)
def ansatz(buffer, t0):
    X(0)
    Ry(t0, 1)
    CNOT(1, 0)

# Run VQE with given ansatz kernel
initAngle = .50

ansatz(buffer, initAngle)

print(buffer.listExtraInfoKeys())
print('Energy = ', buffer.getInformation('vqe-energy'))
print('Opt Angles = ', buffer.getInformation('vqe-angles'))

# Print all children names
print (buffer.getChildrenNames())

# Get all children with the given name
children = buffer.getChildren('Z0')
for c in children:
    # Print any info the children may have 
    print(c.listExtraInfoKeys(), c.getInformation('kernel'), c.getInformation('parameters'))

# Get all children that have 'parameters' information equal to [0.5]
cs = buffer.getChildren('parameters',[0.5])
for c in cs:
    # Print the kernel name and the parameter ([0.5]) and the exp value
    print(c.getInformation('kernel'), c.getInformation('parameters'), c.getExpectationValueZ())
    
# Get all unique children parameters 
ps=buffer.getAllUnique('parameters')
print(ps,'\n', len(ps))

# Get all unique kernel names (should be 5 of them)
print(buffer.getAllUnique('kernel'))

print(buffer)

xacc.Finalize()
