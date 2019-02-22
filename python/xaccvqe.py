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

def getObservableEnergies(buffer, readout=False):
    energies = []
    if readout:
        readout_energies = []
    ps = buffer.getAllUnique('parameters')
    for p in ps:
        energy = 0.0
        readout_energy = 0.0
        for c in buffer.getChildren('parameters', p):
            coeff = c.getInformation('coefficient')
            exp = c.getInformation("exp-val-z")
            e = coeff * exp
            energy += e
            if readout:
                readout_exp = c.getInformation('ro-fixed-exp-val-z')
                ro_e = coeff * readout_exp
                readout_energy += ro_e
        energies.append(energy)
        if readout:
            readout_energies.append(readout_energy)
    if readout:
        return energies, readout_energies
    else:
        return energies

def generateCSV(buffer, file_name, readout=False):
    ps = buffer.getAllUnique('parameters')
    f = open(file_name+".csv", 'w')
    exp_columns = [c.getInformation('kernel') for c in buffer.getChildren('parameters',ps[0])] + ['<E>']
    f.write(str(exp_columns).replace('[','').replace(']','') + '\n')
    for p in ps:
        energy = 0.0
        for c in buffer.getChildren('parameters', p):
            exp = c.getInformation('ro-fixed-exp-val-z') if readout else c.getInformation('exp-val-z')
            energy += exp * c.getInformation('coefficient')
            f.write(str(exp)+',')
        f.write(str(energy)+'\n')
    f.close()

def main(argv=None):
    return

if __name__ == "__main__":
    sys.exit(main())
