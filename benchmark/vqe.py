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
        self.vqe_options_dict = {}
        self.n_qubits = 0
        self.op = None
        
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

        print(xaccOp)
        self.op = xaccOp
        self.n_qubits = n_qubits
        buffer = qpu.createBuffer('q', n_qubits)
        buffer.addExtraInfo('hamiltonian', str(xaccOp))
        buffer.addExtraInfo('ansatz-qasm', ansatz.toString('q').replace('\\n', '\\\\n'))
        pycompiler = xacc.getCompiler('xacc-py')
        buffer.addExtraInfo('ansatz-qasm-py', '\n'.join(pycompiler.translate('q',ansatz).split('\n')[1:]))
            
        if 'n-execs' in inputParams:
            xacc.setOption('sampler-n-execs', inputParams['n-execs'])
            qpu = xacc.getAcceleratorDecorator('improved-sampling', qpu)

        if 'restart-from-file' in inputParams:
            xacc.setOption('vqe-restart-file', inputParams['restart-from-file'])
            qpu = xacc.getAcceleratorDecorator('vqe-restart',qpu)
            qpu.initialize()
                        
        if 'readout-error' in inputParams and inputParams['readout-error']:
            qpu = xacc.getAcceleratorDecorator('ro-error',qpu) 
            
        vqe_opts = {'task': 'vqe', 'accelerator': qpu, 'ansatz':ansatz}
        
        xacc.setOptions(inputParams)

        self.vqe_options_dict = vqe_opts

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
                    fileName = ".persisted_buffer_%s" % (inputParams['accelerator'])
                    file = open(fileName+'.ab', 'w')
                    file.write(str(buffer))
                    file.close()
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
                exp = c.getInformation('ro-fixed-exp-val-z') if c.hasExtraInfoKey('ro-fixed-exp-val-z') else c.getInformation('exp-val-z')
                energy += exp * c.getInformation('coefficient')
                f.write(','+str(exp))
            f.write(','+str(energy)+'\n')
        f.close()

        if 'richardson-extrapolation' in inputParams and inputParams['richardson-extrapolation']:
            angles = buffer.getInformation('vqe-angles')
            qpu = self.vqe_options_dict['accelerator']
            self.vqe_options_dict['accelerator'] = xacc.getAcceleratorDecorator('rich-extrap',qpu)
            self.vqe_options_dict['task'] = 'compute-energy'
            xaccOp = self.op
            self.vqe_options_dict['vqe-params'] = ','.join([str(x) for x in angles])

            fileNames = {r:"%s_%s_%s_%s" % (os.path.splitext(buffer.getInformation('file-name'))[0],
                            buffer.getInformation('accelerator'),
                            'rich_extrap_'+str(r), timestr)+'.csv' for r in [1,3,5,7]}
                            
            nRE_Execs = 2 if not 'rich-extrap-iter' in inputParams else int(inputParams['rich-extrap-iter'])
            if nRE_Execs < 2:
                print('Richardson Extrapolation needs more than 1 execution. Setting to 2.')
                nRE_execs = 2
                
            for r in [1,3,5,7]:
                f = open(fileNames[r], 'w')
                xacc.setOption('rich-extrap-r',r)

                for i in range(nRE_Execs):
                    richardson_buffer = qpu.createBuffer('q', self.n_qubits)
                    results = xaccvqe.execute(xaccOp, richardson_buffer, **self.vqe_options_dict)
                
                    ps = richardson_buffer.getAllUnique('parameters')
                    for p in ps:
                        f.write(str(p).replace('[', '').replace(']', ''))
                        energy = 0.0
                        for c in richardson_buffer.getChildren('parameters', p):
                            exp = c.getInformation('ro-fixed-exp-val-z') if c.hasExtraInfoKey('ro-fixed-exp-val-z') else c.getInformation('exp-val-z')
                            energy += exp * c.getInformation('coefficient')
                            f.write(','+str(exp))
                        f.write(','+str(energy)+'\n')
                f.close()

            nParams = len(ps[0])
            columns = ['t{}'.format(i) for i in range(nParams)]

            kernelNames = [c.getInformation('kernel') for c in buffer.getChildren('parameters',ps[0])]
            columns += kernelNames
            columns.append('E')

            dat = [np.genfromtxt(fileNames[1], delimiter=',', names=columns),
                np.genfromtxt(fileNames[3], delimiter=',', names=columns),
                np.genfromtxt(fileNames[5], delimiter=',', names=columns),
                np.genfromtxt(fileNames[7], delimiter=',', names=columns)]

            allExps = [{k:[] for k in kernelNames} for i in range(4)]
            allEnergies = []

            temp = {r:[] for r in range(4)}

            for i in range(nRE_Execs):
                for r in range(4):
                    for term in kernelNames:
                        allExps[r][term].append(dat[r][term][i])
                    temp[r].append(dat[r]['E'][i])

            evars = [np.std(temp[r]) for r in range(4)]
            xVals = [1,3,5,7]

            avgExps = {k:[np.mean(allExps[r][k]) for r in range(4)] for k in kernelNames}
            varExps = {k:[np.std(allExps[r][k]) for r in range(4)] for k in kernelNames}
            energies = [np.mean(temp[r]) for r in range(4)]

            def f(x, a, b):
                return a*x + b
            res = curve_fit(f, xVals, energies, [1.,energies[0]], sigma=evars)
            print('\nnoisy energy: ', energies[0])
            print('\nrich_extrap intercept: ', res[0][1],'+- ', np.sqrt(np.diag(res[1])[1]))

        if 'hf-energy' in inputParams:
            hf_energy = ast.literal_eval(inputParams['hf-energy'])
            energy = buffer.getInformation('vqe-energy')
            correlation_energy = energy - hf_energy
            buffer.addExtraInfo('correlation-energy', correlation_energy)
