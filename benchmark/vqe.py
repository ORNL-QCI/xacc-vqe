#!usr/bin/python3

from pelix.ipopo.decorators import ComponentFactory, Property, Requires, BindField, UnbindField, Provides, \
    Validate, Invalidate, Instantiate

import ast, configparser, xacc, xaccvqe
from algorithm import Algorithm

@ComponentFactory("vqe_algorithm_factory")
@Provides("xacc_algorithm_service")
@Property("_algorithm","algorithm","vqe")
@Property("_name", "name", "vqe")
@Requires("_molecule_generator", "molecule_generator_service", aggregate=True)
@Instantiate("vqe_algorithm_instance")
class VQE(Algorithm):

    def __init__(self):
        self.molecule_generators = {}
        print("VQE Algorithm...")

    @Validate
    def validate(self, context):
        print("VQE Algorithm Validated")

    @BindField('_molecule_generator')
    def bind_dict(self, field, service, svc_ref):
        print('binding molecule gen')
        generator = svc_ref.get_property('molecule_generator')
        self.molecule_generators[generator] = service

    @UnbindField('_molecule_generator')
    def unbind_dict(self, field, service, svc_ref):

        generator = svc_ref.get_property('molecule_generator')

        del self.molecule_generators[generator]


    @Invalidate
    def invalidate(self, context):
        print("VQE Algorithm Stopped")

    def execute(self, inputParams):
        print("VQE Algorithm Executing")
        qpu = xacc.getAccelerator(inputParams['accelerator'])
        ir_generator = xacc.getIRGenerator(inputParams['ansatz'])
        n_electrons = int(inputParams['n-electrons'])
        xaccOp = xaccvqe.compile(
            self.molecule_generators[inputParams['molecule-generator']].generate(inputParams))
        n_qubits = xaccOp.nQubits()
        buffer = qpu.createBuffer(inputParams['qubit-register'], n_qubits)
        function = ir_generator.generate([n_electrons, n_qubits])
        results = xaccvqe.execute(xaccOp, buffer, **{'task': 'vqe',
                                         'n-electrons': n_electrons})
        return buffer

    def analyze(self, buffer, inputParams):
        print(buffer.getInformation())
