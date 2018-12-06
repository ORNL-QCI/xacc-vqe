from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    BindField, UnbindField, Provides, Validate,
                                    Invalidate, Instantiate)
import ast
import configparser
import xacc
import xaccvqe
from scipy.optimize import minimize
import numpy as np
import time
import os
from xacc import Algorithm
from vqe_base import VQEBase

@ComponentFactory("vqe_algorithm_factory")
@Provides("xacc_algorithm_service")
@Property("_algorithm", "algorithm", "vqe")
@Property("_name", "name", "vqe")
@Requires("_ansatz_generator", "ansatz_generator_service", aggregate=True)
@Requires("_molecule_generator", "molecule_generator_service", aggregate=True)
@Instantiate("vqe_algorithm_instance")
class VQE(VQEBase):

    def __init__(self):
        super().__init__()
        
    @Validate
    def validate(self, context):
        print("VQE Algorithm Validated")

    @Invalidate
    def invalidate(self, context):
        print("VQE Algorithm Invalidated")

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
        
        Adds options:
        - 'scipy-[METHOD]': uses scipy.optimize instead of default nelder-mead to optimize the parameters
    """
        super().execute(inputParams)  
        self.vqe_options_dict['task'] = 'vqe'
        if (inputParams['optimizer'] != 'nelder-mead'):
            if 'scipy' in inputParams['optimizer']:
                scipy_opts = {}
                scipy_opts['method'] = inputParams['optimizer'].replace('scipy-', '')
                if 'tol' in inputParams:
                    scipy_opts['tol'] = ast.literal_eval(inputParams['tol'])
                if 'options' in inputParams:
                    scipy_opts['options'] = ast.literal_eval(inputParams['options'])
                else:
                    scipy_opts['options'] = {'disp': True}

                energies = []
                paramStrings = []
                def energy(p):
                    paramStr = ','.join([str(x) for x in p])
                    e = xaccvqe.execute(self.op, self.buffer, **{'task': 'compute-energy',
                                                           'ansatz': self.ansatz,
                                                           'vqe-params':paramStr,
                                                           'accelerator': self.qpu}).energy
                    energies.append(e)
                    paramStrings.append(paramStr)
                    fileName = ".persisted_buffer_%s" % (inputParams['accelerator'])
                    file = open(fileName+'.ab', 'w')
                    file.write(str(self.buffer))
                    file.close()
                    return e
                if 'initial-parameters' in inputParams:
                    init_params = ast.literal_eval(inputParams['initial-parameters'])
                else:
                    init_params = np.random.uniform(
                        low=-np.pi, high=np.pi, size=(self.ansatz.nParameters(),))
                opt_result = minimize(
                    energy, init_params, **scipy_opts)
                self.buffer.addExtraInfo('vqe-energies',energies)
                self.buffer.addExtraInfo('vqe-parameters',paramStrings)
                return self.buffer
            else:
                xacc.setOption('vqe-backend', inputParams['optimizer'])

        result = xaccvqe.execute(self.op, self.buffer, **self.vqe_options_dict)
        return self.buffer

    def analyze(self, buffer, inputParams):
        super().analyze(buffer, inputParams)