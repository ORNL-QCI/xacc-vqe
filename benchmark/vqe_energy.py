from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    BindField, UnbindField, Provides, Validate,
                                    Invalidate, Instantiate)

import ast
import configparser
import xacc, numpy as np
import xaccvqe
import time
import os
from xacc import Algorithm
from vqe_base import VQEBase
from scipy.optimize import curve_fit

@ComponentFactory("VQE_energy_algorithm_factory")
@Provides("xacc_algorithm_service")
@Property("_algorithm", "algorithm", "vqe-energy")
@Property("_name", "name", "vqe-energy")
@Requires("_hamiltonian_generator", "hamiltonian_generator_service", aggregate=True)
@Requires("_ansatz_generator", "ansatz_generator_service", aggregate=True)
@Instantiate("VQE_energy_algorithm_instance")
class VQEEnergy(VQEBase):
    """
        Algorithm class inherited from VQEBase to execute the 'compute-energy' task of VQE
    """
    def __init__(self):
        super().__init__()

    @Validate
    def validate(self, context):
        """
            iPOPO method that is called when all of the class dependencies have been injected and the class is registered to the framework
        """
        print("VQEEnergy Algorithm Validated")

    @Invalidate
    def invalidate(self, context):
        """
            iPOPO method that is called when the class is removed from the framework or one of its dependencies has been removed
        """
        print("VQEEnergy Algorithm Invalidated")

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
        """
        super().execute(inputParams)
        self.vqe_options_dict['task'] = 'compute-energy'
        results = xaccvqe.execute(self.op, self.buffer, **self.vqe_options_dict)
        return self.buffer

    def analyze(self, buffer, inputParams):
        super().analyze(buffer, inputParams)
