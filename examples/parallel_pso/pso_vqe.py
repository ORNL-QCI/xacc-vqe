import xaccvqe as vqe
import xacc
import numpy as np
from stochopy import Evolutionary
from mpi4py import MPI

# Initialize MPI and XACC
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
xacc.Initialize()
xacc.setOption('n-electrons', '2')

op = vqe.compile("""&FCI NORB=  2,NELEC=  2,MS2= 0,
  ORBSYM=1,5,
  ISYM=0,
 /                         i   a   j   b
  0.6744931033260081E+00   1   1   1   1
  0.6634720448605567E+00   2   2   1   1
  0.6973979494693358E+00   2   2   2   2
  0.1812875358123322E+00   2   1   2   1
 -0.1252477303982147E+01   1   1   0   0
 -0.4759344611440753E+00   2   2   0   0
  0.7137758743754461E+00   0   0   0   0
""")

# Define PSO objective function
def f(params):
    # note 'no-mpi' at the compute energy, hamiltonian term 
    # level, but we are using mpi at the higher PSO swarm level
    return vqe.execute(op,**{'mpi-provider':'no-mpi', 
        'task':'compute-energy', 
        'vqe-params':str(params[0])+','+str(params[1])}).energy

# Construct Stochopy PSO algorithm and execute
ea = Evolutionary(f, popsize=30, 
        lower=np.array([-np.pi,-np.pi]), 
        upper=np.array([np.pi, np.pi]), 
        n_dim=2, mpi=True)
ea.optimize(solver="pso")

# Print result
if mpi_rank == 0:
    print('Our result = ', ea)
