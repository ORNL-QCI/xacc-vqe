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
              getParams = lambda params: ','.join(map(str, params)) 

              execParams = {'ansatz':compiledKernel.getIRFunction(), 'task':'vqe'}
              obs = self.kwargs['observable']
              ars = list(args)
              if len(ars) > 0:
                 arStr = getParams(ars)
                 execParams['vqe-params'] = arStr

              if 'optimizer' in self.kwargs:
                  optimizer = self.kwargs['optimizer']
                  if 'scipy-' in optimizer:
                      optimizer = optimizer.replace('scipy-','')
                      from scipy.optimize import minimize
                      execParams['task'] = 'compute-energy'
                      def energy(params):
                          pStr = getParams(params)
                          execParams['vqe-params'] = pStr
                          return execute(obs, **execParams).energy
                      if len(ars) == 0:
                          import random
                          pi = 3.141592653
                          ars = [random.uniform(-pi,pi) for _ in range(compiledKernel.getIRFunction().nParameters())]
                      opt_result = minimize(energy, ars, method=optimizer, options={'disp':True})
                      return VQETaskResult(opt_result.fun, opt_result.x)
                  else:
                      xacc.setOption('vqe-backend', optimizer)
                      
              return execute(obs, **execParams)
          return wrapped_f

def main(argv=None):
   return
   
if __name__ == "__main__":
    sys.exit(main())
