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


@ComponentFactory("vqe_algorithm_factory")
@Provides("xacc_algorithm_service")
@Property("_algorithm", "algorithm", "vqe")
@Property("_name", "name", "vqe")
@Requires("_ansatz_generator", "ansatz_generator_service", aggregate=True)
@Requires("_molecule_generator", "molecule_generator_service", aggregate=True)
@Instantiate("vqe_algorithm_instance")
class VQE(Algorithm):

    def __init__(self):
        # define the list of MoleculeGenerator services installed and available
        self.molecule_generators = {}
        self.ansatz_generators = {}

    @Validate
    def validate(self, context):
        print("VQE Algorithm Validated")

    @Invalidate
    def invalidate(self, context):
        print("VQE Algorithm Invalidated")

    @BindField('_ansatz_generator')
    @BindField('_molecule_generator')
    def bind_dicts(self, field, service, svc_ref):
        if svc_ref.get_property('molecule_generator'):
            generator = svc_ref.get_property('molecule_generator')
            self.molecule_generators[generator] = service
        elif svc_ref.get_property('ansatz_generator'):
            generator = svc_ref.get_property('ansatz_generator')
            self.ansatz_generators[generator] = service

    @UnbindField('_ansatz_generator')
    @UnbindField('_molecule_generator')
    def unbind_dicts(self, field, service, svc_ref):

        if svc_ref.get_property('molecule_generator'):
            generator = svc_ref.get_property('molecule_generator')
            del self.molecule_generators[generator]
        elif svc_ref.get_property('ansatz_generator'):
            generator = svc_ref.get_property('ansatz_generator')
            del self.ansatz_generators[generator]

    def execute(self, inputParams):
        qpu = xacc.getAccelerator(inputParams['accelerator'])
        xaccOp = self.molecule_generators[inputParams['molecule-generator']].generate(
            inputParams)
        ansatz = self.ansatz_generators[inputParams['name']].generate(
            inputParams, xaccOp.nQubits())

        if 'qubit-map' in inputParams:
            qubit_map = ast.literal_eval(inputParams['qubit-map'])
            xaccOp, ansatz, n_qubits = xaccvqe.mapToPhysicalQubits(
                xaccOp, ansatz, qubit_map)
        else:
            n_qubits = xaccOp.nQubits()

        buffer = qpu.createBuffer('q', n_qubits)
        buffer.addExtraInfo('hamiltonian', str(xaccOp))
 
        if 'readout-error' in inputParams and inputParams['readout-error']:
            qpu = xacc.getAcceleratorDecorator('ro-error',qpu)
  
        vqe_opts = {'task': 'vqe', 'accelerator': qpu, 'ansatz':ansatz}

        xacc.setOptions(inputParams)

        if 'initial-parameters' in inputParams:
            vqe_opts['vqe-params'] = inputParams['initial-parameters']
        
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
                    e = xaccvqe.execute(xaccOp, buffer, **{'task': 'compute-energy',
                                                           'ansatz': ansatz,
                                                           'vqe-params':paramStr,
                                                           'accelerator': qpu}).energy
                    energies.append(e)
                    paramStrings.append(paramStr)
                    return e
                if 'initial-parameters' in inputParams:
                    init_params = ast.literal_eval(inputParams['initial-parameters'])
                else:
                    init_params = np.random.uniform(
                        low=-np.pi, high=np.pi, size=(ansatz.nParameters(),))
                print(str(scipy_opts))
                opt_result = minimize(
                    energy, init_params, **scipy_opts)
                buffer.addExtraInfo('vqe-energies',energies)
                buffer.addExtraInfo('vqe-parameters',paramStrings)
                return buffer
            else:
                xacc.setOption('vqe-backend', inputParams['optimizer'])

        result = xaccvqe.execute(xaccOp, buffer, **vqe_opts)
        return buffer

    def analyze(self, buffer, inputParams):
        ps = buffer.getAllUnique('parameters')
        timestr = time.strftime("%Y%m%d-%H%M%S")
        csv_name = "%s_%s_%s_%s" % (os.path.splitext(buffer.getInformation('file-name'))[0],
                                    buffer.getInformation('accelerator'),
                                    'vqe', timestr)
        f = open(csv_name+".csv", 'w')
        for p in ps:
            f.write(str(p).replace('[', '').replace(']', ''))
            energy = 0.0
            for c in buffer.getChildren('parameters', p):
                exp = c.getInformation('exp-val-z')
                energy += exp * c.getInformation('coefficient')
                f.write(','+str(c.getInformation('exp-val-z')))
            f.write(','+str(energy)+'\n')
        f.close()
