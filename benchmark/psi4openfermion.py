from pelix.ipopo.decorators import ComponentFactory, Property, Requires, Provides, \
    Validate, Invalidate, Instantiate
from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4
from openfermion.transforms import get_fermion_operator
from moleculegenerator import MoleculeGenerator
import ast


@ComponentFactory("Psi4OF_molecule_generator_factory")
@Provides("molecule_generator_service")
@Property("_molecule_generator", "molecule_generator", "psi4Of")
@Property("_name", "name", "psi4Of")
@Instantiate("Psi4OF_molecule_generator_instance")
class Psi4OpenFermion(MoleculeGenerator):

    def __init__(self):
        print("MoleculeGenerator initialized")

    @Validate
    def validate(self, context):
        print("MoleculeGenerator validated")

    @Invalidate
    def invalidate(self, context):
        print("MoleculeGenerator invalidated")

    def generate(self, inputParams):
        mdata = MolecularData(ast.literal_eval(inputParams['geometry']),
                              inputParams['basis'],
                              int(inputParams['multiplicity']),
                              int(inputParams['charge']))
        molecule = run_psi4(mdata, run_scf=1, run_fci=1)
        fermiOp = get_fermion_operator(molecule.get_molecular_hamiltonian())
        return fermiOp
