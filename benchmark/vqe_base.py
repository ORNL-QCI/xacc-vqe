from xacc import Algorithm
import xacc
import xaccvqe
import ast
import time
import os
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