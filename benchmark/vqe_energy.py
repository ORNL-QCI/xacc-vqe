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
@Requires("_molecule_generator", "molecule_generator_service", aggregate=True)
@Requires("_ansatz_generator", "ansatz_generator_service", aggregate=True)
@Instantiate("VQE_energy_algorithm_instance")
class VQEEnergy(VQEBase):

    def __init__(self):
        super().__init__()
        
    @Validate
    def validate(self, context):
        print("VQEEnergy Algorithm Validated")

    @Invalidate
    def invalidate(self, context):
        print("VQEEnergy Algorithm Invalidated")
    
    @BindField('_ansatz_generator')
    @BindField('_molecule_generator')
    def bind_dicts(self, field, service, svc_ref):
        super().bind_dicts(field, service, svc_ref)

    @UnbindField('_ansatz_generator')
    @UnbindField('_molecule_generator')
    def unbind_dicts(self, field, service, svc_ref):
        super().unbind_dicts(field, service, svc_ref)

    def execute(self, inputParams):
    """
        Inherited method with algorithm-specific implementation
        
        - sets XACC VQE task to 'compute-energy' 
    """
        super().execute(inputParams)
        self.vqe_options_dict['task'] = 'compute-energy'
        results = xaccvqe.execute(self.op, self.buffer, **self.vqe_options_dict)
        return self.buffer

    def analyze(self, buffer, inputParams):
        super().analyze(buffer, inputParams) 