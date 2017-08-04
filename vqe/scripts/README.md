# How to use this script

To run this script, you need Psi4 and Fermilib (which also gets you ProjectQ). For Psi4, 

```bash
$ git clone https://github.com/psi4/psi4
$ cd psi4 && mkdir build && cd build
$ cmake ..
$ make install
$ export PYTHONPATH=/usr/local/psi4/lib:$PYTHONPATH
```

For Fermilib,

```bash
$ pip install fermilib --global-option=--without-cppsimulator
```

Then you should be able to run this script as follows

```bash
$ python generateVQEKernels -m test_molecules/h2.mol -a .392
```

This will generate a *.hpp file that contains the XACC Kernel representation 
of this molecule's fermionic hamiltonian. 
