#ifndef VQETASKS_COMPUTEENERGYVQETASK_HPP_
#define VQETASKS_COMPUTEENERGYVQETASK_HPP_

#include "VQETask.hpp"

namespace xacc {
namespace vqe {

class ComputeEnergyVQETask : public VQETask {

public:
  ComputeEnergyVQETask() {}

  ComputeEnergyVQETask(std::shared_ptr<VQEProgram> prog) : VQETask(prog) {}

  virtual VQETaskResult execute(Eigen::VectorXd parameters);

  /**
   * Return the name of this instance.
   *
   * @return name The string name
   */
  virtual const std::string name() const { return "compute-energy"; }

  /**
   * Return the description of this instance
   * @return description The description of this object.
   */
  virtual const std::string description() const {
    return "This VQETask computes the energy at the given set of parameters.";
  }

  /**
   * Return an empty options_description, this is for
   * subclasses to implement.
   */
  virtual OptionPairs getOptions() {
    OptionPairs desc {{"vqe-use-mpi", "Use MPI distributed execution."},{
        "vqe-persist-data",
        "Base file name for buffer data."},{
        "converge-ro-error", "Use ro-fixed-exp-val-z to compute energy."}};
    return desc;
  }

  int vqeIteration = 0;
  int totalQpuCalls = 0;
};
} // namespace vqe
} // namespace xacc
#endif
