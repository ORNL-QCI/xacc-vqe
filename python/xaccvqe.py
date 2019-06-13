from _pyxaccvqe import *
import xacc
import os
import platform
import sys
import sysconfig
import argparse
import inspect
from abc import abstractmethod, ABC
from xacc import PauliOperator

class VQEOpt(ABC):

    @abstractmethod
    def optimize(self, observable, buffer, optimizer_args, execParams):
        self.opt_args = optimizer_args
        self.execParams = execParams
        self.energies = []
        self.obs = observable
        self.buffer = buffer

    @abstractmethod
    def energy(self, params):
        pStr = ",".join(map(str, params))
        self.execParams['vqe-params'] = pStr
        e = execute(self.obs, self.buffer, **self.execParams).energy
        self.energies.append(e)
        return e

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

def getPurifiedEnergies(buffer, raw=False):
    ps = buffer.getAllUnique('parameters')
    p_es = []
    nonp_es = []
    for p in ps:
        for c in buffer.getChildren('parameters', p):
            if c.name() == 'I':
                continue
            pure_energy = c.getInformation('purified-energy')
            if raw:
                nonpure_energy = c.getInformation('non-purified-energy')
        p_es.append(pure_energy)
        if raw:
            nonp_es.append(nonpure_energy)
    if raw:
        return p_es, nonp_es
    else:
        return p_es

def variance(singleChildBuffer, nPhysicalBits):
    import numpy as np
    c = 0.0
    op = PauliOperator()
    name = singleChildBuffer.name()
    # FIXME, 2 is hardcoded, should be computed dynamically
    bufAsString = " ".join(name[i:i+2] for i in range(0, len(name), 2))
    op.fromString('(1,0) '+bufAsString)

    data = singleChildBuffer.getMeasurementCounts()
    exp = singleChildBuffer.getExpectationValueZ()
    z,x = op.toBinaryVectors(nPhysicalBits)
    pzorx = np.logical_or(z,x)

    nshots = 0
    for k, v in data.items():
        bitstr = np.asarray(list(k))[::-1].astype(np.bool)
        s = -1.0 if np.logical_xor.reduce(np.logical_and(bitstr, pzorx)) else 1.0
        c += (s - exp ) * (s - exp ) * v
        nshots += v
    c /= (nshots - 1)
    return c

def covariance(singleChildBufferA, singleChildBufferB, nPhysicalBits):
    import numpy as np
    c = 0.0
    opA = PauliOperator()
    opB = PauliOperator()
    nameA = singleChildBufferA.name()
    nameB = singleChildBufferB.name()

    bufAsStringA = " ".join(nameA[i:i+2] for i in range(0, len(nameA), 2))
    opA.fromString('(1,0) '+bufAsStringA)
    bufAsStringB = " ".join(nameB[i:i+2] for i in range(0, len(nameB), 2))
    opB.fromString('(1,0) '+bufAsStringB)

    data = singleChildBufferA.getMeasurementCounts()
    exp = singleChildBufferA.getExpectationValueZ()
    za,xa = opA.toBinaryVectors(nPhysicalBits)
    pzorxa = np.logical_or(za,xa)
    zb,xb = opB.toBinaryVectors(nPhysicalBits)
    pzorxb = np.logical_or(zb,xb)

    nshots = 0
    for k, v in data.items():
        bitstr = np.asarray(list(k))[::-1].astype(np.bool)
        sa = -1.0 if np.logical_xor.reduce(np.logical_and(bitstr, pzorxa)) else 1.0
        sb = -1.0 if np.logical_xor.reduce(np.logical_and(bitstr, pzorxb)) else 1.0
        c += (sa - exp ) * (sb - exp ) * v
        nshots += v
    c /= (nshots - 1)
    return c

def main(argv=None):
    return

if __name__ == "__main__":
    sys.exit(main())
