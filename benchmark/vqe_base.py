from xacc import Algorithm
import xacc
import xaccvqe
import ast
import time
import os
import numpy as np 

class VQEBase(Algorithm):
    
    def __init__(self):
        self.molecule_generators = {}
        self.ansatz_generators = {}
        self.vqe_options_dict = {}
        self.n_qubits = 0
        self.buffer = None
        self.ansatz = None
        self.op = None
        self.qpu = None
        
    def execute(self, inputParams):
        self.qpu = xacc.getAccelerator(inputParams['accelerator'])
        xaccOp = self.molecule_generators[inputParams['molecule-generator']].generate(
            inputParams)
        self.ansatz = self.ansatz_generators[inputParams['name']].generate(
            inputParams, xaccOp.nQubits())
        if 'qubit-map' in inputParams:
            qubit_map = ast.literal_eval(inputParams['qubit-map'])
            xaccOp, self.ansatz, n_qubits = xaccvqe.mapToPhysicalQubits(
                xaccOp, self.ansatz, qubit_map)
        else:
            n_qubits = xaccOp.nQubits()
        self.op = xaccOp
        self.n_qubits = n_qubits
        self.buffer = self.qpu.createBuffer('q', n_qubits)
        self.buffer.addExtraInfo('hamiltonian', str(xaccOp))
        self.buffer.addExtraInfo('ansatz-qasm', self.ansatz.toString('q').replace('\\n', '\\\\n'))
        pycompiler = xacc.getCompiler('xacc-py')
        #self.buffer.addExtraInfo('ansatz-qasm-py', '\n'.join(pycompiler.translate('q',self.ansatz).split('\n')[1:]))
        if 'n-execs' in inputParams:
            xacc.setOption('sampler-n-execs', inputParams['n-execs'])
            self.qpu = xacc.getAcceleratorDecorator('improved-sampling', self.qpu)

        if 'restart-from-file' in inputParams:
            xacc.setOption('vqe-restart-file', inputParams['restart-from-file'])
            self.qpu = xacc.getAcceleratorDecorator('vqe-restart',self.qpu)
            self.qpu.initialize()
                        
        if 'readout-error' in inputParams and inputParams['readout-error']:
            self.qpu = xacc.getAcceleratorDecorator('ro-error',self.qpu) 
            
        self.vqe_options_dict = {'accelerator': self.qpu, 'ansatz': self.ansatz}
        
        if 'initial-parameters' in inputParams:
            self.vqe_options_dict['vqe-params'] = ','.join([str(x) for x in ast.literal_eval(inputParams['initial-parameters'])])
            
        xacc.setOptions(inputParams)            
    
    def analyze(self, buffer, inputParams):
        ps = buffer.getAllUnique('parameters')
        timestr = time.strftime("%Y%m%d-%H%M%S")
        exp_csv_name = "%s_%s_%s_%s" % (os.path.splitext(buffer.getInformation('file-name'))[0],
                                        buffer.getInformation('accelerator'),"exp_val_z",
                                        timestr)
        f = open(exp_csv_name+".csv", 'w')
        exp_columns = [c.getInformation('kernel') for c in buffer.getChildren('parameters',ps[0])] + ['<E>']
        f.write(str(exp_columns).replace('[','').replace(']','') + '\n')
        for p in ps:
            energy = 0.0
            for c in buffer.getChildren('parameters', p):
                exp = c.getInformation('exp-val-z')
                energy += exp * c.getInformation('coefficient')
                f.write(str(exp)+',')
            f.write(str(energy)+'\n')
        f.close()
        ## Repeating code - putting in one loop would be messy, but possible
        if 'readout-error' in inputParams:
            ro_exp_csv_name = "%s_%s_%s_%s" % (os.path.splitext(buffer.getInformation('file-name'))[0],
                                        buffer.getInformation('accelerator'),"ro_fixed_exp_val_z",
                                        timestr)
            f = open(ro_exp_csv_name+'.csv', 'w')
            f.write(str(exp_columns).replace('[','').replace(']','')+'\n')
            for p in ps:
                energy = 0.0
                for c in buffer.getChildren('parameters', p):
                    exp = c.getInformation('ro-fixed-exp-val-z')
                    energy += exp * c.getInformation('coefficient')
                    f.write(str(exp)+',')
                f.write(str(energy)+'\n')
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

            def linear(x, a, b):
                return a*x+b
    
            def exp(x, a,b):
                return a*np.exp(b*x)# + b

            from scipy.optimize import curve_fit
            res = curve_fit(linear, xVals, energies, [1,energies[0]], sigma=evars)

            print(res)
            print('\nnoisy energy: ', energies[0])
            print('\nrich linear extrap: ', res[0][1],'+- ', np.sqrt(np.diag(res[1])[1]))

            res_exp = curve_fit(exp, xVals, energies, [0,0], sigma=evars)

            print('\nrich exp extrap: ', exp(0,res_exp[0][0],res_exp[0][1]), '+-', np.sqrt(np.diag(res_exp[1])[1]))

            print('\n')
            print('E(r): ', energies)
            print('a*exp(b*r):', [exp(x,res_exp[0][0],res_exp[0][1]) for x in xVals])
