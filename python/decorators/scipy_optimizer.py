from pelix.ipopo.decorators import (ComponentFactory, Property, Instantiate)
import xacc
from xaccvqe import VQEOpt
import xaccvqe as vqe
import inspect
from abc import abstractmethod, ABC
from scipy.optimize import minimize

@ComponentFactory("scipy_opt_factory")
@Property("_vqe_optimizer", "vqe_optimizer", "scipy-opt")
@Instantiate("scipy_opt_instance")
class ScipyOpt(VQEOpt):

    def optimize(self, observable, buffer, optimizer_args, execParams):
        super().optimize(observable, buffer, optimizer_args, execParams)

        opt_result = minimize(self.energy, self.init_args, **self.opt_args)

        # Optimizer adds the results to the buffer automatically
        buffer.addExtraInfo('vqe-energies', self.energies)
        buffer.addExtraInfo('vqe-parameters', self.angles)
        optimal_angles = [float(x) for x in self.angles[self.energies.index(min(self.energies))].split(",")]
        buffer.addExtraInfo('vqe-angles', optimal_angles)
        buffer.addExtraInfo('vqe-energy', min(self.energies))

    # Noticing something very weird if the objective energy function
    # resides in the super class; looks like it needs to be
    # redefined everytime (in an optimizer)
    def energy(self, params):
        pStr = ",".join(map(str, params))
        self.execParams['vqe-params'] = pStr
        e = vqe.execute(self.obs, self.buffer, **self.execParams).energy

        if 'rdm-purification' in self.execParams['accelerator'].name():
            t = self.buffer.getAllUnique('parameters')
            ind = len(t) - 1
            children = self.buffer.getChildren('parameters', t[ind])
            e = children[1].getInformation('purified-energy')

        self.angles.append(self.execParams['vqe-params'])
        self.energies.append(e)
        fileName = ".persisted_buffer_%s" % (self.buffer.getInformation('accelerator')) if self.buffer.hasExtraInfoKey('accelerator') \
                                                                else ".persisted_buffer_%s" % (self.execParams['accelerator'].name())
        file = open(fileName+'.ab', 'w')
        file.write(str(self.buffer))
        file.close()
        return e



