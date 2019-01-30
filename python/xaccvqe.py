from _pyxaccvqe import *
import xacc
import os
import platform
import sys
import sysconfig
import argparse
import inspect

def QubitOperator2XACC(qubit_op):
    """
    Map openfermion qubit operator to an XACC PauliOperator
    """
    xaccOp = PauliOperator()
    for k, v in qubit_op.terms.items():
        terms = dict(k)
        if terms:
            xaccOp += PauliOperator(terms, v)
        else:
            xaccOp += PauliOperator(v)
    return(xaccOp)

def XACC2QubitOperator(pauli_op):
    try:
        from openfermion import QubitOperator
    except:
        xacc.error("OpenFermion not installed, cannot convert PauliOperator to QubitOperator.")
        return

    qop = QubitOperator()
    for o in pauli_op:
        term = tuple(o[1].ops().items())
        if term == () or term[0][1] == 'I':
            qop += QubitOperator((), o[1].coeff())
        else:
            qop += QubitOperator(tuple(o[1].ops().items()), o[1].coeff())
    return qop

def mapToPhysicalQubits(op, ansatz, logical2PhysicalMap):
    n_qubits = max(logical2PhysicalMap) + 1
    ir = op.toXACCIR()
    xacc.setOption('qubit-map',','.join([str(i) for i in logical2PhysicalMap]))
    irp = xacc.getIRPreprocessor('qubit-map-preprocessor')
    irp.process(ir)

    ham = PauliOperator()
    ham.fromXACCIR(ir)

    ansatzir = xacc.gate.createIR()
    ansatzir.addKernel(ansatz)
    irp.process(ansatzir)
    xacc.unsetOption('qubit-map')
    return ham, ansatz, n_qubits

def main(argv=None):
    return

if __name__ == "__main__":
    sys.exit(main())
