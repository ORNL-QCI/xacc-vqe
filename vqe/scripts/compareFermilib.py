import scipy
import scipy.linalg
from scipy.optimize import minimize
import random
import glob, os, sys
import numpy as np

# fermilib/projectQ
from fermilib.ops import FermionOperator
from fermilib.transforms import get_fermion_operator, get_sparse_operator, jordan_wigner
from fermilib.utils import get_ground_state, MolecularData, uccsd_trotter_engine, uccsd_singlet_evolution, uccsd_singlet_paramsize, uccsd_singlet_operator
from fermilibpluginpsi4 import run_psi4
from projectq.ops import X, All, Measure
from projectq.backends import CommandPrinter, CircuitDrawer
import psi4

def H2(r):
    """ di-hydrogen molecule with seperation length r

    Parameters
    ----------
    r: interatomic distance in units of Angstroms

    Retruns
    -------
    mol: psi4.core.Molecule type representation of H2(r)
    """

    mol = psi4.geometry("""
    0 1
    H
    H 1 {0}
    """.format(r)
    )
    return mol

def fl_geo(mol):
    """
    output molecule geometry for fermilib cartesian representation
    Parameters:
    -----------
    mol: psi4.core.Molecule object

    Returns:
    --------
    geometry tuple with lower case second letters
    and XYZ coordinates for each atom
    """
    mol.update_geometry()
    np_geo = np.array(mol.geometry())
    return [(mol.label(i)[0] + str.lower(str(mol.label(i)[1:])), tuple(np_geo[i])) for i in range(mol.natom())]

def f(packed_amplitudes):
   compiler_engine = uccsd_trotter_engine()
   wavefunction = compiler_engine.allocate_qureg(4)
   for i in range(2):
      X | wavefunction[i]
   evolution_op = uccsd_singlet_evolution(packed_amplitudes, 4, 2)
   evolution_op | wavefunction
   compiler_engine.flush()
   energy = compiler_engine.backend.get_expectation_value(qubit_hamiltonian, wavefunction)
   print "expectation = ", energy
   All(Measure) | wavefunction
   compiler_engine.flush()
   return energy

def printCircuit(packed_amplitudes):
   compiler_engine = uccsd_trotter_engine(CommandPrinter())
   wavefunction = compiler_engine.allocate_qureg(4)
   for i in range(2):
      X | wavefunction[i]
   evolution_op = uccsd_singlet_evolution(packed_amplitudes, 4, 2)
   evolution_op | wavefunction
   compiler_engine.flush()
   return

diatomic_bond_length = 1.4 #1645885
newmol = H2(diatomic_bond_length)
geometry = fl_geo(newmol)
basis = 'sto-3g'
multiplicity = newmol.multiplicity()
charge = newmol.molecular_charge()
description = str(diatomic_bond_length)
molecule = MolecularData(geometry, basis, multiplicity,
                         charge, description)
molecule.save()

molecule = run_psi4(molecule,run_scf=1,run_mp2=1,
                    run_cisd=0,run_ccsd=0,run_fci=1)
molecular_hamiltonian = molecule.get_molecular_hamiltonian()
fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
fermion_hamiltonian.compress()
print('The fermionic Hamiltonian in canonical basis follows:\n{}'.format(fermion_hamiltonian))
qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
qubit_hamiltonian.compress()
print('The Jordan-Wigner Hamiltonian in canonical basis follows:\n{}'.format(qubit_hamiltonian))

initial_amplitudes = [0, .05677]
op = jordan_wigner(uccsd_singlet_operator(initial_amplitudes, 4, 2))
print "UCCSD:\n", op
printCircuit(initial_amplitudes)
print f(initial_amplitudes)

opt_result = minimize(f, initial_amplitudes,
                      method="CG", options={'disp':True})

opt_energy, opt_amplitudes = opt_result.fun, opt_result.x
print "Bond Distance: ", diatomic_bond_length
print("\nOptimal UCCSD Singlet Energy: {}".format(opt_energy))
print("Optimal UCCSD Singlet Amplitudes: {}".format(opt_amplitudes))
print("Classical CCSD Energy: {} Hartrees".format(molecule.ccsd_energy))
print("Exact FCI Energy: {} Hartrees".format(molecule.fci_energy))
