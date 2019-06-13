from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    Provides, Instantiate)
import xacc
from xaccvqe import VQEOpt
import xaccvqe as vqe
import inspect
from abc import abstractmethod, ABC
import pybobyqa
import numpy as np

@ComponentFactory("bobyqa_opt_factory")
@Provides("vqe_optimization")
@Property("_name", "name", "bobyqa-opt")
@Property("_vqe_optimizer", "vqe_optimizer", "bobyqa-opt")
@Instantiate("bobyqa_opt_instance")
class BOBYQAOpt(VQEOpt):

    def optimize(self, observable, buffer, optimizer_args, execParams):
        super().optimize(observable, buffer, optimizer_args, execParams)

        if 'vqe-params' in self.execParams:
            init_args = [float(x) for x in self.execParams['vqe-params'].split(',')]
        else:
            import random
            pi = 3.141592653
            init_args = np.array([random.uniform(-pi, pi) for _ in range(self.execParams['ansatz'].nParameters())])

        opt_result = pybobyqa.solve(self.energy, init_args, **self.opt_args)

    # For some reason, the Py-BOBYQA module
    # will not work if using super().energy(params) ( like in ScipyOpt )
    # so, had to redefine here.
    def energy(self, params):
        pStr = ",".join(map(str, params))
        self.execParams['vqe-params'] = pStr
        e = vqe.execute(self.obs, self.buffer, **self.execParams).energy
        self.angles.append(self.execParams['vqe-params'])
        self.energies.append(e)
        return e

