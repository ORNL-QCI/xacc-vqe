from pelix.ipopo.decorators import ComponentFactory, Property, Requires, Provides, \
    Validate, Invalidate, Instantiate
from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4
from openfermion.transforms import get_fermion_operator
from moleculegenerator import MoleculeGenerator
import ast
from xaccvqe import PauliOperator


@ComponentFactory("customH2_molecule_generator_factory")
@Provides("molecule_generator_service")
@Property("_molecule_generator", "molecule_generator", "customH2")
@Property("_name", "name", "customH2")
@Instantiate("customH2_molecule_generator_instance")
class CustomH2(MoleculeGenerator):

    def __init__(self):
        print("MoleculeGenerator initialized")

    @Validate
    def validate(self, context):
        print("MoleculeGenerator validated")

    @Invalidate
    def invalidate(self, context):
        print("MoleculeGenerator invalidated")

    def generate(self, inputParams):
        return PauliOperator(.2976) + PauliOperator({0:'Z'},.3593) + PauliOperator({1:'Z'},-.4826) + PauliOperator({0:'Z',1:'Z'},.5818) + PauliOperator({0:'Y',1:'Y'}, .0896) + PauliOperator({0:'X',1:'X'},.0896)
