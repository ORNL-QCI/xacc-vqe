from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    Provides, Instantiate)
import xacc
from xaccvqe import VQEOpt
import xaccvqe as vqe
import inspect
from abc import abstractmethod, ABC
from scipy.optimize import minimize

@ComponentFactory("scipy_opt_factory")
@Provides("vqe_optimization")
@Property("_name", "name", "scipy-opt")
@Property("_vqe_optimizer", "vqe_optimizer", "scipy-opt")
@Instantiate("scipy_opt_instance")
class ScipyOpt(VQEOpt):

    def optimize(self, observable, buffer, optimizer_args, execParams):
        super().optimize(observable, buffer, optimizer_args, execParams)
        print(optimizer_args)
        print(self.opt_args)
        self.opt_args['options'] = {"maxiter": 30}
        if 'vqe-params' in self.execParams:
            init_args = [float(x) for x in self.execParams['vqe-params'].split(',')]
        else:
            import random
            pi = 3.141592653
            init_args = [random.uniform(-pi, pi) for _ in range(self.execParams['ansatz'].nParameters())]

        opt_result = minimize(self.energy, init_args, **self.opt_args)

    def energy(self, params):
        super().energy(params)


