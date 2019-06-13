from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    Provides, Instantiate, BindField, UnbindField)
import xacc
import inspect
import xaccvqe as vqe
from _pyxaccvqe import *

@ComponentFactory("wrapped_vqe_factory")
@Provides("decorator_algorithm_service")
@Property("_algorithm", "algorithm", "vqe")
@Property("_name", "name", "vqe")
@Requires("_vqe_optimizers", "vqe_optimization", aggregate=True)
@Instantiate("wrapped_vqe_instance")
class WrappedVQEF(xacc.DecoratorFunction):

    def __init__(self):
        self.vqe_optimizers = {}

    @BindField("_vqe_optimizers")
    def bind_optimizers(self, field, service, svc_ref):
        if svc_ref.get_property('vqe_optimizer'):
            optimizer = svc_ref.get_property('vqe_optimizer')
            self.vqe_optimizers[optimizer] = service

    @UnbindField("_vqe_optimizers")
    def unbind_optimizers(self, field, service, svc_ref):
        if svc_ref.get_property('vqe_optimizer'):
            optimizer = svc_ref.get_property('vqe_optimizer')
            del vqe_optimizers[optimizer]

    def __call__(self, *args, **kwargs):
        super().__call__(*args, **kwargs)
        def getParams(params): return ','.join(map(str, params))
        execParams = {'accelerator': self.qpu, 'ansatz': self.compiledKernel, 'task': 'vqe'}
        obs = self.kwargs['observable']
        ars = list(args)

        if not isinstance(args[0], xacc.AcceleratorBuffer):
            raise RuntimeError(
                'First argument of an xacc kernel must be the Accelerator Buffer to operate on.')

        buffer = ars[0]
        ars = ars[1:]
        if len(ars) > 0:
            arStr = getParams(ars)
            execParams['vqe-params'] = arStr
        if 'optimizer' in self.kwargs:
            opt_name = self.kwargs['optimizer']
            optimizer_args = {}
            optimizer = self.vqe_optimizers[opt_name]
            if 'options' in self.kwargs:
                optimizer_args = self.kwargs['options']
            optimizer.optimize(obs, buffer, optimizer_args, execParams)
            buffer.addExtraInfo('vqe-energies', optimizer.energies)
            # we want to use scipy
            # if 'scipy-' in optimizer:
            #     execParams['task'] = 'compute-energy'
            #     # get the scipy-opt optimizer plugin
            #     opt_instance = self.vqe_optimizers['scipy-opt']
            #     optimizer = optimizer.replace('scipy-', '')
            #     # set the scipy-optimizer to the specified method
            #     optimizer_args["method"] = optimizer
            #     if 'options' in self.kwargs:
            #         optimizer_args['options'] = self.kwargs['options']
            #     # call the optimizer-plugin optimize() method
            #     # with the optimizer settings and VQE settings
            #     opt_instance.optimize(obs, buffer, optimizer_args, execParams)

                # now we don't need any of this

                # from scipy.optimize import minimize
                # energies = []
                # def energy(params):
                #     pStr = getParams(params)
                #     execParams['vqe-params'] = pStr
                #     e = vqe.execute(obs, buffer, **execParams).energy
                #     energies.append(e)
                #     return e
                # if len(ars) == 0:
                #     import random
                #     pi = 3.141592653
                #     ars = [
                #         random.uniform(-pi, pi) for _ in range(self.compiledKernel.nParameters())]
                # optargs = {'method': optimizer, 'options': {'disp': True}}
                # if 'options' in self.kwargs:
                #     optargs['options'] = self.kwargs['options']
                # if 'tol' in optargs['options']:
                #     print(optargs['options']['tol'])
                # if 'opt_params' in self.kwargs:
                #     for k, v in self.kwargs['opt_params'].items():
                #         optargs[k] = v

                # opt_result = minimize(energy, ars, **optargs)
                # print(optargs)
                # print(opt_result)

                # the optimizer has all of the energies now
                # buffer.addExtraInfo('vqe-energies',opt_instance.energies)
                # return
            # until we have a standard naming for these too
            # use bobyqa instead

            # elif 'bobyqa' or 'BOBYQA' in optimizer:
            #     execParams['task'] = 'compute-energy'
            #     opt_instance = self.vqe_optimizers['bobyqa-opt']
            #     if 'options' in self.kwargs:
            #         optimizer_args = self.kwargs['options']
            #     opt_instance.optimize(obs, buffer, optimizer_args, execParams)

                # we dont need any of this anymore

                # import pybobyqa
                # energies = []
                # def energy(params):
                #     pStr = getParams(params)
                #     execParams['vqe-params'] = pStr
                #     e = vqe.execute(obs, buffer, **execParams).energy
                #     return e
                # if len(ars) == 0:
                #     import random
                #     pi = 3.141592653
                #     ars = np.array([
                #         random.uniform(-pi, pi) for _ in range(self.compiledKernel.nParameters())])
                # print(type(ars))
                # print(ars)
                # opt_result = pybobyqa.solve(energy, ars)
                # print(opt_result)

                # buffer.addExtraInfo('vqe-energies', opt_instance.energies)
        #     else:
        #         xacc.setOption('vqe-backend', optimizer)
        #         if 'opt_params' in self.kwargs:
        #             for k, v in self.kwargs['opt_params'].items():
        #                 xacc.setOption(k, str(v))
        # vqe.execute(obs, buffer, **execParams)
        return