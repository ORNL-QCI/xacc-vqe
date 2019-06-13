from pelix.ipopo.decorators import (ComponentFactory, Property, Requires,
                                    Provides, Instantiate, BindField, UnbindField)
import xacc
import inspect
import xaccvqe as vqe

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
        # optimizer given
        if 'optimizer' in self.kwargs:
            opt_name = self.kwargs['optimizer']
            optimizer_args = {}
            # get optimizer from available vqe_optimization services
            optimizer = self.vqe_optimizers[opt_name]
            # get all the options to pass to optimizer
            if 'options' in self.kwargs:
                optimizer_args = self.kwargs['options']
            # call optimize() method of optimizer
            optimizer.optimize(obs, buffer, optimizer_args, execParams)
            # add the energies obtained and the angles with which they were obtained
            # to the AcceleratorBuffer
            # Also, add the optimal angles to the buffer
            buffer.addExtraInfo('vqe-energies', optimizer.energies)
            buffer.addExtraInfo('vqe-parameters', optimizer.angles)
            vqe_angles = [float(x) for x in optimizer.angles[optimizer.energies.index(min(optimizer.energies))].split(",")]
            buffer.addExtraInfo('vqe-angles', vqe_angles)
        else:
            vqe.execute(obs, buffer, **execParams)
        return