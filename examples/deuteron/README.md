## Cloud Quantum Computing of an Atomic Nucleus 

Here we present the code used to compute the binding energy of deuteron on 
both IBM and Rigetti quantum computers. This directory gives a couple examples 
of leveraging XACC and the XACC-VQE application to compute the binding energy of deuteron. 

First, we have provided 2 Jupyter notebooks, one for the H2 Hamiltonian and 
another for the H3 Hamiltonian. Each targets the TNQVM MPS simulator 
by default, but users can change this by specifying the desired Accelerator 
in the notebooks (update the string from 'tnqvm' to 'ibm' for instance). 

We also provide a description of how to run these problems with the 
XACC-VQE command line executable. You can find the kernel source code 
and a description of the command to run in the command_line directory. 

