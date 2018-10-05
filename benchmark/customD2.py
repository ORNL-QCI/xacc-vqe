from pelix.ipopo.decorators import ComponentFactory, Property, Requires, Provides, \
    Validate, Invalidate, Instantiate
from moleculegenerator import MoleculeGenerator
import ast
from xaccvqe import PauliOperator


@ComponentFactory("customD2_molecule_generator_factory")
@Provides("molecule_generator_service")
@Property("_molecule_generator", "molecule_generator", "customD2")
@Property("_name", "name", "customD2")
@Instantiate("customD2_molecule_generator_instance")
class CustomD2(MoleculeGenerator):

    def __init__(self):
        print("MoleculeGenerator initialized")

    @Validate
    def validate(self, context):
        print("MoleculeGenerator validated")

    @Invalidate
    def invalidate(self, context):
        print("MoleculeGenerator invalidated")

    def generate(self, inputParams):
        return PauliOperator(5.906709445) + PauliOperator({0:'X',1:'X'}, -2.1433) + PauliOperator({0:'Y',1:'Y'}, -2.1433) + PauliOperator({0:'Z'}, .21829) + PauliOperator({1:'Z'}, -6.125)
