from pelix.ipopo.decorators import (ComponentFactory, Property, Instantiate)
import ast
import xacc
import xaccvqe
from xacc import BenchmarkAlgorithm
from vqe_base import VQEBase

@ComponentFactory("vqe_benchmark_factory")
@Property("_name", "name", "vqe")
@Instantiate("vqe_benchmark")
class VQE(VQEBase):
    """
        Algorithm class inherited from VQEBase to execute the VQE algorithm
    """
    def __init__(self):
        super().__init__()

    def execute(self, inputParams):
        """
            Inherited method with algorithm-specific implementation

            Parameters:
                inputParams - a dictionary of input parameters obtained from .ini file
            Adds options:
            - 'scipy-[METHOD]': uses scipy.optimize instead of default nelder-mead to optimize the parameters
        """
        super().execute(inputParams)
        self.vqe_options_dict['task'] = 'vqe'
        if self.optimizer is not None:
            self.optimizer.optimize(self.op, self.buffer, self.optimizer_options, self.vqe_options_dict)
        else:
            result = xaccvqe.execute(self.op, self.buffer, **self.vqe_options_dict)

        return self.buffer

    def analyze(self, buffer, inputParams):
        super().analyze(buffer, inputParams)
