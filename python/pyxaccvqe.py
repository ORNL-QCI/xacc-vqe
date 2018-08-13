from _pyxaccvqe import *
import pyxacc as xacc
import os, platform, sys, sysconfig
import argparse, inspect

class qpu(xacc.qpu):
   class vqe(object):
       def __init__(self, *args, **kwargs):
           self.args = args
           self.kwargs = kwargs
           return
       def __call__(self, f):
          def wrapped_f(*args, **kwargs):
              src = '\n'.join(inspect.getsource(f).split('\n')[1:])
              compiler = xacc.getCompiler('xacc-py')
              if isinstance(self.kwargs['accelerator'], xacc.Accelerator):
                  qpu = self.kwargs['accelerator']
              else:
                  qpu = xacc.getAccelerator(self.kwargs['accelerator'])
              ir = compiler.compile(src, qpu)
              program = xacc.Program(qpu, ir)
              compiledKernel = program.getKernels()[0]
              obs = self.kwargs['observable']
              ars = list(args)
              if len(ars) > 0:
                 arStr = str(ars[0])
                 for i in ars[1:]:
                    arStr += ','+str(i)
                 return execute(obs, **{'ansatz':compiledKernel.getIRFunction(), 'vqe-params':arStr, 'task':'vqe'})
              else:
                 return execute(obs, **{'ansatz':compiledKernel.getIRFunction(), 'task':'vqe'})

          return wrapped_f

def main(argv=None):
   return
   
if __name__ == "__main__":
    sys.exit(main())
