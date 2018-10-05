from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    BindField, UnbindField, Provides, Validate, 
                                    Invalidate, Instantiate)

import ast
import configparser
import xacc
import xaccvqe
from algorithm import Algorithm

@ComponentFactory("VQE_energy_algorithm_factory")
@Provides("xacc_algorithm_service")
@Property("_algorithm", "algorithm", "vqe-energy")
@Property("_name", "name", "vqe-energy")
@Requires("_molecule_generator", "molecule_generator_service", aggregate=True)
@Instantiate("VQE_energy_algorithm_instance")
class VQEEnergy(Algorithm):

    def __init__(self):
        
        self.molecule_generators = {}

    @Validate
    def validate(self, context):
        print("VQEEnergy Algorithm Validated")

    @BindField('_molecule_generator')
    def bind_dict(self, field, service, svc_ref):

        generator = svc_ref.get_property('molecule_generator')
        self.molecule_generators[generator] = service

    @UnbindField('_molecule_generator')
    def unbind_dict(self, field, service, svc_ref):

        generator = svc_ref.get_property('molecule_generator')

        del self.molecule_generators[generator]

    @Invalidate
    def invalidate(self, context):
        print("VQEEnergy Algorithm Invalidated")

    def execute(self, inputParams):
        qpu = xacc.getAccelerator(inputParams['accelerator'])
        ir_generator = xacc.getIRGenerator(inputParams['ansatz'])
        xaccOp = self.molecule_generators[inputParams['molecule-generator']].generate(inputParams)
        n_qubits = xaccOp.nQubits()
        buffer = qpu.createBuffer('q', n_qubits)
        function = ir_generator.generate(ast.literal_eval(inputParams['ansatz-params']))
        results = xaccvqe.execute(xaccOp, buffer, **{'task': 'compute-energy',
                                                     'ansatz': function,
                                                     'vqe-params': inputParams['vqe-params']})
        return buffer

    def analyze(self, buffer, inputParams):
        print(buffer.getInformation())
