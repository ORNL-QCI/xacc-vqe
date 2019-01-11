from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    Provides, Instantiate)
import xacc
import inspect
from _pyxaccvqe import *

@ComponentFactory("wrapped_vqe_factory")
@Provides("decorator_algorithm_service")
@Property("_algorithm", "algorithm", "vqe")
@Property("_name", "name", "vqe")
@Instantiate("wrapped_vqe_instance")
class WrappedVQEF(xacc.DecoratorFunction):

    def __call__(self, *args, **kwargs):
        super().__call__(*args, **kwargs)

        def getParams(params): return ','.join(map(str, params))
        execParams = {'accelerator': self.qpu, 'ansatz': self.compiledKernel.getIRFunction(), 'task': 'vqe'}
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
                        random.uniform(-pi, pi) for _ in range(self.compiledKernel.getIRFunction().nParameters())]
                optargs = {'method': optimizer, 'options': {'disp': True}}
                if 'options' in self.kwargs:
                    optargs['options'] = self.kwargs['options']
                if 'opt_params' in self.kwargs:
                    for k, v in self.kwargs['opt_params'].items():
                        optargs[k] = v
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