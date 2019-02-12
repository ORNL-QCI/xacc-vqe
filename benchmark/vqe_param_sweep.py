from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    BindField, UnbindField, Provides, Validate,
                                    Invalidate, Instantiate)

import ast
import configparser
import xacc, numpy as np
import xaccvqe
import time
import os
from xacc import BenchmarkAlgorithm
from vqe_base import VQEBase
from _pyxaccvqe import *

@ComponentFactory("VQE_param_sweep_algorithm_factory")
@Provides("benchmark_algorithm_service")
@Property("_algorithm", "algorithm", "param-sweep")
@Property("_name", "name", "param-sweep")
@Requires("_hamiltonian_generator", "hamiltonian_generator_service", aggregate=True)
@Requires("_ansatz_generator", "ansatz_generator_service", aggregate=True)
@Instantiate("VQE_param_sweep_algorithm_instance")
class ParamSweep(VQEBase):
    """
        Algorithm class inherited from VQEBase to execute the 'compute-energy' task of VQE
    """
    def __init__(self):
        super().__init__()

    @BindField('_ansatz_generator')
    @BindField('_hamiltonian_generator')
    def bind_dicts(self, field, service, svc_ref):
        """
            iPOPO method to bind ansatz and hamiltonian generator dependencies for use by the Algorithm bundle
        """
        super().bind_dicts(field, service, svc_ref)

    @UnbindField('_ansatz_generator')
    @UnbindField('_hamiltonian_generator')
    def unbind_dicts(self, field, service, svc_ref):
        """
            iPOPO method to unbind ansatz and hamiltonian generator dependencies for use by the Algorithm bundle

            Called when the bundle is invalidated
        """
        super().unbind_dicts(field, service, svc_ref)

    def execute(self, inputParams):
        """
            Inherited method with algorithm-specific implementation

            Parameters:
                inputParams - a dictionary of input parameters obtained from .ini file

            - sets XACC VQE task to 'compute-energy'
            - executes a parameter-sweep
        """
        super().execute(inputParams)
        self.vqe_options_dict['task'] = 'compute-energy'
        
        if inputParams['upper-bound'] == 'pi':
            up_bound = np.pi
        else:
            up_bound = ast.literal_eval(inputParams['upper-bound'])

        if inputParams['lower-bound'] == '-pi':
            low_bound = -np.pi
        else:
            low_bound = ast.literal_eval(inputParams['lower-bound']) 
        
        num_params = ast.literal_eval(inputParams['num-params'])
        
        for param in np.linspace(low_bound, up_bound, num_params):
            self.vqe_options_dict['vqe-params'] = str(param)
            results = execute(self.op, self.buffer, **self.vqe_options_dict)
        return self.buffer

    def analyze(self, buffer, inputParams):
        super().analyze(buffer, inputParams)
