from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    BindField, UnbindField, Provides, Validate,
                                    Invalidate, Instantiate)
import xacc
import inspect
from _pyxaccvqe import *

@ComponentFactory("wrapped_energy_factory")
@Provides("decorator_algorithm_service")
@Property("_algorithm", "algorithm", "energy")
@Property("_name", "name", "energy")
@Instantiate("wrapped_energy_instance")
class WrappedEnergyF(xacc.DecoratorFunction):

    def __call__(self, *args, **kwargs):
        super().__call__(*args, **kwargs)

        def getParams(params): return ','.join(map(str, params))

        execParams = {'accelerator': self.qpu, 'ansatz': self.compiledKernel.getIRFunction(), 'task': 'compute-energy'}
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