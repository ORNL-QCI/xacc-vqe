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

def mapToPhysicalQubits(op, ansatz, logical2PhysicalMap):
    n_qubits = max(logical2PhysicalMap) + 1
    ir = op.toXACCIR()
    xacc.setOption('qubit-map',','.join([str(i) for i in [3,4]]))
    irp = xacc.getIRPreprocessor('qubit-map-preprocessor')
    irp.process(ir)

    ham = PauliOperator()
    ham.fromXACCIR(ir)
    
    ansatzir = xacc.gate.createIR()
    ansatzir.addKernel(ansatz)
    irp.process(ansatzir)
    xacc.unsetOption('qubit-map')
    return ham, ansatz, n_qubits
    
class WrappedVQEF(xacc.WrappedF):
    def __init__(self, f, *args, **kwargs):
        xacc.WrappedF.__init__(self, f, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        src = '\n'.join(inspect.getsource(self.function).split('\n')[1:])
        compiler = xacc.getCompiler('xacc-py')
        if isinstance(self.kwargs['accelerator'], xacc.Accelerator):
            qpu = self.kwargs['accelerator']
        else:
            qpu = xacc.getAccelerator(self.kwargs['accelerator'])
        ir = compiler.compile(src, qpu)
        program = xacc.Program(qpu, ir)
        compiledKernel = program.getKernels()[0]

        def getParams(params): return ','.join(map(str, params))
        execParams = {'accelerator': qpu, 'ansatz': compiledKernel.getIRFunction(), 'task': 'vqe'}
        obs = self.kwargs['observable']
        ars = list(args)

        if not isinstance(args[0], xacc.AcceleratorBuffer):
            raise RuntimeError(
                'First argument of an xacc kernel must be the Accelerator Buffer to operate on.')

        buffer = ars[0]
        ars = ars[1:]
        if len(ars) > 0:
            arStr = getParams(ars)
            execParams['vqe-params'] = arStr

        if 'optimizer' in self.kwargs:
            optimizer = self.kwargs['optimizer']
            if 'scipy-' in optimizer:
                optimizer = optimizer.replace('scipy-', '')
                from scipy.optimize import minimize
                execParams['task'] = 'compute-energy'

                energies = []
                def energy(params):
                    pStr = getParams(params)
                    execParams['vqe-params'] = pStr
                    e = execute(obs, buffer, **execParams).energy
                    energies.append(e)
                    return e
                if len(ars) == 0:
                    import random
                    pi = 3.141592653
                    ars = [
                        random.uniform(-pi, pi) for _ in range(compiledKernel.getIRFunction().nParameters())]
                optargs = {'method': optimizer, 'options': {'disp': True}}
                if 'options' in self.kwargs:
                    print(self.kwargs['options'])
                    optargs['options'] = self.kwargs['options']
                if 'opt_params' in self.kwargs:
                    for k, v in self.kwargs['opt_params'].items():
                        optargs[k] = v
                print(optargs)
                opt_result = minimize(energy, ars, **optargs)
                buffer.addExtraInfo('vqe-energies',energies)
                return
            else:
                xacc.setOption('vqe-backend', optimizer)
                if 'opt_params' in self.kwargs:
                    for k, v in self.kwargs['opt_params'].items():
                        xacc.setOption(k, str(v))
        execute(obs, buffer, **execParams)
        return


class WrappedEnergyF(xacc.WrappedF):
    def __init__(self, f, *args, **kwargs):
        xacc.WrappedF.__init__(self, f, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        src = '\n'.join(inspect.getsource(self.function).split('\n')[1:])
        compiler = xacc.getCompiler('xacc-py')
        if isinstance(self.kwargs['accelerator'], xacc.Accelerator):
            qpu = self.kwargs['accelerator']
        else:
            qpu = xacc.getAccelerator(self.kwargs['accelerator'])
        ir = compiler.compile(src, qpu)
        program = xacc.Program(qpu, ir)
        compiledKernel = program.getKernels()[0]

        def getParams(params): return ','.join(map(str, params))

        execParams = {'accelerator': qpu, 'ansatz': compiledKernel.getIRFunction(), 'task': 'compute-energy'}
        obs = self.kwargs['observable']
        ars = list(args)

        if not isinstance(args[0], xacc.AcceleratorBuffer):
            raise RuntimeError(
                'First argument of an xacc kernel must be the Accelerator Buffer to operate on.')

        buffer = ars[0]
        ars = ars[1:]
        if len(ars) > 0:
            arStr = getParams(ars)
            execParams['vqe-params'] = arStr
        execute(obs, buffer, **execParams)
        return


class qpu(xacc.qpu):
    class vqe(object):
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
            self.__dict__.update(kwargs)
            return

        def __call__(self, f):
            return WrappedVQEF(f, *self.args, **self.kwargs)

    class energy(object):
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
            self.__dict__.update(kwargs)
            return

        def __call__(self, f):
            return WrappedEnergyF(f, *self.args, **self.kwargs)


def main(argv=None):
    return


if __name__ == "__main__":
    sys.exit(main())
