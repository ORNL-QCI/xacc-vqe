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
        super().execute(inputParams)
        self.vqe_options_dict['task'] = 'compute-energy'
        results = xaccvqe.execute(self.op, self.buffer, **self.vqe_options_dict)
        return self.buffer

    def analyze(self, buffer, inputParams):
        super().analyze(buffer, inputParams)
#         ps = buffer.getAllUnique('parameters')
#         timestr = time.strftime("%Y%m%d-%H%M%S")
#         csv_name = "%s_%s_%s_%s" % (os.path.splitext(buffer.getInformation('file-name'))[0],
#                                         buffer.getInformation('accelerator'),
#                                         'vqeenergy', timestr)
#         f = open(csv_name+".csv", 'w')
#         for p in ps:
#             f.write(str(p).replace('[', '').replace(']', ''))
#             energy = 0.0
#             for c in buffer.getChildren('parameters', p):
#                 exp = c.getInformation('exp-val-z')
#                 energy += exp * c.getInformation('coefficient')
#                 f.write(','+str(c.getInformation('exp-val-z')))
#             f.write(','+str(energy)+'\n')
#         f.close()
       
#         if 'richardson-extrapolation' in inputParams and inputParams['richardson-extrapolation']:
#             angles = buffer.getInformation('vqe-angles')
#             qpu = self.vqe_options_dict['accelerator']
#             self.vqe_options_dict['accelerator'] = xacc.getAcceleratorDecorator('rich-extrap',qpu)
#             self.vqe_options_dict['task'] = 'compute-energy'
#             xaccOp = self.op
#             self.vqe_options_dict['vqe-params'] = ','.join([str(x) for x in angles])

#             fileNames = {r:"%s_%s_%s_%s" % (os.path.splitext(buffer.getInformation('file-name'))[0],
#                             buffer.getInformation('accelerator'),
#                             'rich_extrap_'+str(r), timestr)+'.csv' for r in [1,3,5,7]}
                            
#             nRE_Execs = 2 if not 'rich-extrap-iter' in inputParams else int(inputParams['rich-extrap-iter'])
#             if nRE_Execs < 2:
#                 print('Richardson Extrapolation needs more than 1 execution. Setting to 2.')
#                 nRE_execs = 2
                
#             for r in [1,3,5,7]:
#                 f = open(fileNames[r], 'w')
#                 xacc.setOption('rich-extrap-r',r)

#                 for i in range(nRE_Execs):
#                     richardson_buffer = qpu.createBuffer('q', self.n_qubits)
#                     results = xaccvqe.execute(xaccOp, richardson_buffer, **self.vqe_options_dict)
                
#                     ps = richardson_buffer.getAllUnique('parameters')
#                     for p in ps:
#                         f.write(str(p).replace('[', '').replace(']', ''))
#                         energy = 0.0
#                         for c in richardson_buffer.getChildren('parameters', p):
#                             exp = c.getInformation('ro-fixed-exp-val-z') if c.hasExtraInfoKey('ro-fixed-exp-val-z') else c.getInformation('exp-val-z')
#                             energy += exp * c.getInformation('coefficient')
#                             f.write(','+str(exp))
#                         f.write(','+str(energy)+'\n')
#                 f.close()

#             nParams = len(ps[0])
#             columns = ['t{}'.format(i) for i in range(nParams)]

#             kernelNames = [c.getInformation('kernel') for c in buffer.getChildren('parameters',ps[0])]
#             columns += kernelNames
#             columns.append('E')

#             dat = [np.genfromtxt(fileNames[1], delimiter=',', names=columns),
#                 np.genfromtxt(fileNames[3], delimiter=',', names=columns),
#                 np.genfromtxt(fileNames[5], delimiter=',', names=columns),
#                 np.genfromtxt(fileNames[7], delimiter=',', names=columns)]

#             allExps = [{k:[] for k in kernelNames} for i in range(4)]
#             allEnergies = []

#             temp = {r:[] for r in range(4)}

#             for i in range(nRE_Execs):
#                 for r in range(4):
#                     for term in kernelNames:
#                         allExps[r][term].append(dat[r][term][i])
#                     temp[r].append(dat[r]['E'][i])

#             evars = [np.std(temp[r]) for r in range(4)]
#             xVals = [1,3,5,7]

#             avgExps = {k:[np.mean(allExps[r][k]) for r in range(4)] for k in kernelNames}
#             varExps = {k:[np.std(allExps[r][k]) for r in range(4)] for k in kernelNames}
#             energies = [np.mean(temp[r]) for r in range(4)]

#             def f(x, a, b):
#                 return a*x + b
#             res = curve_fit(f, xVals, energies, [1.,energies[0]], sigma=evars)
#             print('\nnoisy energy: ', energies[0])
#             print('\nrich_extrap intercept: ', res[0][1],'+- ', np.sqrt(np.diag(res[1])[1]))
            
#         if 'hf-energy' in inputParams:
#             hf_energy = ast.literal_eval(inputParams['hf-energy'])
#             energy = buffer.getInformation('vqe-energy')
#             correlation_energy = energy - hf_energy
#             buffer.addExtraInfo('correlation-energy', correlation_energy)        