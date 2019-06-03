from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    BindField, UnbindField, Provides, Validate,
                                    Invalidate, Instantiate)
import ast
import configparser
import xacc
import xaccvqe
import time
import os
from xacc import BenchmarkAlgorithm
from vqe_base import VQEBase

@ComponentFactory("param-sweep_benchmark_factory")
@Provides("benchmark_algorithm")
@Property("_name", "name", "param-sweep")
@Requires("_hamiltonian_generators", "hamiltonian_generator", aggregate=True)
@Requires("_ansatz_generators", "ansatz_generator", aggregate=True)
@Instantiate("param-sweep_benchmark")
class ParamSweep(VQEBase):
    """
        Algorithm class inherited from VQEBase to execute the 'compute-energy' task of VQE
    """
    def __init__(self):
        super().__init__()

    @BindField('_ansatz_generators')
    @BindField('_hamiltonian_generators')
    def bind_dicts(self, field, service, svc_ref):
        """
            iPOPO method to bind ansatz and hamiltonian generator dependencies for use by the Algorithm bundle
        """
        super().bind_dicts(field, service, svc_ref)

    @UnbindField('_ansatz_generators')
    @UnbindField('_hamiltonian_generators')
    def unbind_dicts(self, field, service, svc_ref):
        """
            iPOPO method to unbind ansatz and hamiltonian generator dependencies for use by the Algorithm bundle

            Called when the bundle is invalidated
        """
        super().unbind_dicts(field, service, svc_ref)

    def linspace(self, a, b, n=100):
        if n < 2:
            return b
        diff = (float(b) - a)/(n - 1)
        return [diff * i + a  for i in range(n)]

    def execute(self, inputParams):
        """
            Inherited method with algorithm-specific implementation

            Parameters:
                inputParams - a dictionary of input parameters obtained from .ini file

            - sets XACC VQE task to 'compute-energy'
            - executes a parameter-sweep
        """
        super().execute(inputParams)
        pi = 3.141592653589793
        self.vqe_options_dict['task'] = 'compute-energy'

        if inputParams['upper-bound'] == 'pi':
            up_bound = pi
        else:
            up_bound = ast.literal_eval(inputParams['upper-bound'])

        if inputParams['lower-bound'] == '-pi':
            low_bound = -pi
        else:
            low_bound = ast.literal_eval(inputParams['lower-bound'])

        num_params = ast.literal_eval(inputParams['num-params'])

        for param in self.linspace(low_bound, up_bound, num_params):
            self.vqe_options_dict['vqe-params'] = str(param)
            results = xaccvqe.execute(self.op, self.buffer, **self.vqe_options_dict)
        return self.buffer

    def analyze(self, buffer, inputParams):
        super().analyze(buffer, inputParams)
