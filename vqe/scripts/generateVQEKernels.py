import numpy as np
import argparse
import psi4
import sys
import os 

from fermilib.ops import FermionOperator
from fermilib.utils import MolecularData, uccsd_singlet_operator
from fermilibpluginpsi4 import run_psi4
from fermilib.transforms import get_fermion_operator, jordan_wigner
from projectq.backends import CommandPrinter, CircuitDrawer
from scipy.optimize import minimize

#psi4.set_memory('2.5 GB')

#psi4.set_options({'reference': 'uhf'})
#psi4.set_options({'scf_type': 'pk'})
#psi4.set_options({'basis': 'sto-3g'})

def parse_args(args):
    """ Parse command line arguments and return them. """
    parser = argparse.ArgumentParser(description="XACC VQE Fermion Kernel Generator.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     fromfile_prefix_chars='@')
    parser.add_argument("-m", "--molecule", required=True)
    parser.add_argument("-a", "--molecule-args", nargs='*', type=float, help="The arguments for the molecule generation source string.")
    parser.add_argument("-r", "--args-range", nargs='*', type=float, help="The start, end, and step for a range of args")

    opts = parser.parse_args(args)
    return opts

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
    return [(mol.label(i)[0] + str.lower(str(mol.label(i)[1:]))
             , tuple(np_geo[i])) for i in range(mol.natom())]

def main(argv=None):
    opts = parse_args(sys.argv[1:])

    moleculeFile = opts.molecule
    args = opts.molecule_args
    r = opts.args_range
	
    print r
	
    if (args == None and r != None):
        args = list(np.arange(r[0],r[1],r[2]))

    src = ''
    with open(moleculeFile, 'r') as myfile:
       src=myfile.read()

    print args
    print src
    exec(src)

    for arg in args:
       print 'generating arg: ', arg
       molecule = generateMolecule(arg)

       geometry = fl_geo(molecule)
       print 'Geom: ', geometry
       basis = 'sto-3g'
       multiplicity = molecule.multiplicity()
       charge = molecule.molecular_charge()
       description = str('H2 Molecule')

       # Make molecule and print out a few interesting facts about it.
       moleculeData = MolecularData(geometry, basis, multiplicity,
                            charge, description)
       moleculeData.save()

       print('Molecule has automatically generated name {}'.format(
           moleculeData.name))
       print('Information about this molecule would be saved at:\n{}\n'.format(
           moleculeData.filename))
       print('This molecule has {} atoms and {} electrons.'.format(
           moleculeData.n_atoms, moleculeData.n_electrons))

       for atom, atomic_number in zip(moleculeData.atoms, moleculeData.protons):
           print('Contains {} atom, which has {} protons.'.format(
               atom, atomic_number))

       mol = run_psi4(moleculeData,run_scf=1,run_mp2=1,
                           run_cisd=0,run_ccsd=0,run_fci=1)
                        
       # Load full molecular Hamiltonian 
       molecular_hamiltonian = mol.get_molecular_hamiltonian()
       print molecular_hamiltonian
       # Map operator to fermions and qubits.
       fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
       fermion_hamiltonian.compress()
       print fermion_hamiltonian    
       xaccKernelStr = '__qpu__ ' + moleculeData.name.replace(" ","_") + str(arg).replace(".","_") + '() {\n'
       for i, term in enumerate(list(fermion_hamiltonian.terms.keys())): 
          xaccKernelStr += '\t' + str(fermion_hamiltonian.terms[term]) + ' '
          for j, op in enumerate(term):
              xaccKernelStr += str(op[0]) + ' ' + str(op[1]) + ' '
          xaccKernelStr += '\n'

       xaccKernelStr += '}'
       print 'Kernel\n', xaccKernelStr
    
       dirname = moleculeData.name.replace(" ","_")
       if not os.path.exists(dirname):
          os.makedirs(dirname)
       filename = dirname + '_' + str(arg)
       kernelFile = open(dirname+'/'+filename+'.hpp', "w")
       kernelFile.write(xaccKernelStr)
       kernelFile.close()
    
       qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
       qubit_hamiltonian.compress()
       print('The Jordan-Wigner Hamiltonian in canonical basis follows:\n{}'.format(qubit_hamiltonian))
    
       singOp = uccsd_singlet_operator([1,1],4,2)
       print('The UCCSD Singlet Operator follows:\n{}'.format(singOp))
       jwSingOp = jordan_wigner(singOp)
       print('The UCCSD Singlet Operator JW follows:\n{}'.format(jwSingOp))

if __name__ == "__main__":
    sys.exit(main())
