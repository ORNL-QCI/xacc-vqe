#ifndef XACC_PROGRAM_HPP_
#define XACC_PROGRAM_HPP_

#include "Compiler.hpp"
#include "Accelerator.hpp"

namespace xacc {
/**
 * Utility class used to map variadic arguments to
 * a vector of InstructionParameters.
 */
class GetInstructionParametersFunctor {
protected:
  std::vector<InstructionParameter> &params;

public:
  GetInstructionParametersFunctor(std::vector<InstructionParameter> &p)
      : params(p) {}
  template <typename T> void operator()(const T &t) {
    params.push_back(InstructionParameter(t));
  }
};

/**
 * The Kernel represents a functor to be executed
 * on the attached Accelerator.
 */
template <typename... RuntimeArgs> class Kernel {

protected:
  /**
   * The IR Function that this Kernel is executing.
   */
  std::shared_ptr<Function> function;

  /**
   * The Accelerator to execute this Kernel on
   */
  std::shared_ptr<Accelerator> accelerator;

public:
  /**
   * The Constructor
   * @param acc Accelerator to execute on
   * @param f Function to execute
   */
  Kernel(std::shared_ptr<Accelerator> acc, std::shared_ptr<Function> f)
      : function(f), accelerator(acc) {}

  /**
   * The copy constructor
   * @param k Other kernel to copy
   */
  Kernel(const Kernel &k) : function(k.function), accelerator(k.accelerator) {}

  /**
   * Execute this Kernel on the given AcceleratorBuffer and
   * with the provided list of runtime arguments.
   *
   * @param buffer AcceleratorBuffer to execute on.
   * @param args Runtime parameters
   */
  void operator()(std::shared_ptr<AcceleratorBuffer> buffer,
                  RuntimeArgs... args) {
    if (sizeof...(RuntimeArgs) > 0) {
      // Store the runtime parameters in a tuple
      auto argsTuple = std::make_tuple(args...);

      // Loop through the tuple, and add InstructionParameters
      // to the parameters vector.
      std::vector<InstructionParameter> parameters;
      GetInstructionParametersFunctor f(parameters);
      xacc::tuple_for_each(argsTuple, f);

      int count = 0;
      std::vector<double> paramVec(parameters.size());
      for (auto i : parameters) {
        if (i.which() == 0) {
          paramVec[count] = (double)mpark::get<int>(i);
        } else if (i.which() == 1) {
          paramVec[count] = mpark::get<double>(i);
        } else {
          XACCLogger::instance()->error(
              "Kernel.operator() mapping runtime params to Eigen::Vector, "
              "invalid runtime param type.");
        }
        count++;
      }

      // Evaluate all Variable Parameters
      auto evaled = function->operator()(paramVec);
      accelerator->execute(buffer, evaled);
    } else {
      // Execute the Kernel on the Accelerator
      accelerator->execute(buffer, function);
    }
  }

  /**
   * Execute this Kernel on the given AcceleratorBuffer and with
   * the provided vector of InstructionParameters.
   *
   * @param buffer The AcceleratorBuffer
   * @param parameters The runtime parameters
   */
  void operator()(std::shared_ptr<AcceleratorBuffer> buffer,
                  std::vector<InstructionParameter> parameters) {
    int count = 0;
    std::vector<double> paramVec(parameters.size());
    for (auto i : parameters) {
      if (i.which() == 0) {
        paramVec[count] = (double)mpark::get<int>(i);
      } else if (i.which() == 1) {
        paramVec[count] = mpark::get<double>(i);
      } else {
        XACCLogger::instance()->error(
            "Kernel.operator() mapping runtime params to Eigen::Vector, "
            "invalid runtime param type.");
      }
      count++;
    }

    // Evaluate all Variable Parameters
    auto evaled = function->operator()(paramVec);
    accelerator->execute(buffer, evaled);
  }

  /**
   * Return the number of Kernel runtime parameters
   * @return nParams The number of Parameters
   */
  const int getNumberOfKernelParameters() { return function->nParameters(); }

  /**
   * Return the IR Function wrapped by this Kernel.
   *
   * @return function The IR Function
   */
  std::shared_ptr<Function> getIRFunction() { return function; }

  const std::string getName() { return function->name(); }
};

/**
 * The KernelList is a standard C++ vector that provides
 * and overloaded operator() operator that delegates to the
 * Accelerator.execute() method that executes multiple
 * IR Functions in a single execute() call.
 */
template <typename... RuntimeArgs>
class KernelList : public std::vector<Kernel<RuntimeArgs...>> {

protected:
  /**
   * The Accelerator to execute this Kernel on
   */
  std::shared_ptr<Accelerator> accelerator;

public:
  KernelList() {}
  /**
   * The Constructor
   */
  KernelList(std::shared_ptr<Accelerator> acc) : accelerator(acc) {}


  /**
   * Return the Accelerator this KernelList delegates to.
   */
  std::shared_ptr<Accelerator> getAccelerator() { return accelerator; }


  /**
   * Templated operator() overload.
   */
  std::vector<std::shared_ptr<AcceleratorBuffer>>
  operator()(std::shared_ptr<AcceleratorBuffer> buffer, RuntimeArgs... args) {

    std::vector<std::shared_ptr<Function>> functions;
    if (sizeof...(RuntimeArgs) > 0) {
      // Store the runtime parameters in a tuple
      auto argsTuple = std::make_tuple(args...);

      // Loop through the tuple, and add InstructionParameters
      // to the parameters vector.
      std::vector<InstructionParameter> parameters;
      GetInstructionParametersFunctor f(parameters);
      xacc::tuple_for_each(argsTuple, f);

      int count = 0;
      std::vector<double> paramVec(parameters.size());
      for (auto i : parameters) {
        if (i.which() == 0) {
          paramVec[count] = (double)mpark::get<int>(i);
        } else if (i.which() == 1) {
          paramVec[count]= mpark::get<double>(i);
        } else {
          XACCLogger::instance()->error(
              "Kernel.operator() mapping runtime params to Eigen::Vector, "
              "invalid runtime param type.");
        }
        count++;
      }

      for (auto k : *this) {
        auto eval = k.getIRFunction()->operator()(paramVec);
        functions.push_back(eval);
      }
    }

    auto buffers = accelerator->execute(buffer, functions);

    return buffers;
  }

  /**
   * Overloaded operator() operator that takes InstructionParameters.
   */
  std::vector<std::shared_ptr<AcceleratorBuffer>>
  operator()(std::shared_ptr<AcceleratorBuffer> buffer,
             std::vector<InstructionParameter> parameters) {

    int count = 0;
    std::vector<double> paramVec(parameters.size());
    for (auto i : parameters) {
      if (i.which() == 0) {
        paramVec[count] = (double)mpark::get<int>(i);
      } else if (i.which() == 1) {
        paramVec[count] = mpark::get<double>(i);
      } else {
        XACCLogger::instance()->error(
            "Kernel.operator() mapping runtime params to Eigen::Vector, "
            "invalid runtime param type.");
      }
      count++;
    }

    std::vector<std::shared_ptr<Function>> functions;
    for (auto k : *this) {
      auto eval = k.getIRFunction()->operator()(paramVec);
      functions.push_back(eval);
    }

    auto buffers = accelerator->execute(buffer, functions);


    return buffers;
  }

  std::vector<std::shared_ptr<AcceleratorBuffer>>
  execute(std::shared_ptr<AcceleratorBuffer> buffer) {
    std::vector<std::shared_ptr<Function>> functions;
    for (auto k : *this) {
      if (k.getIRFunction()->nInstructions() > 0) {
        functions.push_back(k.getIRFunction());
      }
    }
    auto buffers = accelerator->execute(buffer, functions);

    return buffers;
  }
};

using AcceleratorType = Accelerator::AcceleratorType;

/**
 * Placeholder, soon we will have Acc Independent transformations...
 * @param accType
 * @return
 */
std::vector<IRTransformation>
getAcceleratorIndependentTransformations(AcceleratorType &accType);

/**
 * The Program is the main entrypoint for the XACC
 * API. Users with accelerator kernels must construct a
 * valid Program to be compiled and executed on the
 * attached accelerator. Programs must be given the
 * Accelerator reference to be used and kernel source
 * code at construction time.
 */
class Program {

protected:
  /**
   * Reference to the source accelerator
   * kernel code to be compiled and executed
   */
  std::string src;

  /**
   * Reference to the attached Accelerator to
   * use in this compilation and execution
   */
  std::shared_ptr<Accelerator> accelerator;

  /**
   * Reference to the XACC IR instance that is
   * created by the Compiler
   */
  std::shared_ptr<IR> xaccIR;

  /**
   * Reference to the compiler
   */
  std::shared_ptr<Compiler> compiler;

public:
  /**
   * The Constructor, takes the Accelerator
   * to execute on, and the source to compile and execute
   *
   * @param acc Attached Accelerator to execute
   * @param sourceFile The kernel source code
   */
  Program(std::shared_ptr<Accelerator> acc, const std::string &sourceFile)
      : src(sourceFile), accelerator(std::move(acc)) {
  }

  Program(std::shared_ptr<Accelerator> acc, std::shared_ptr<xacc::IR> ir)
      : xaccIR(ir), accelerator(std::move(acc)) {
  }

  /**
   * The Constructor, takes the Accelerator to execute on,
   * and the source kernel string as a stream.
   * @param acc Attached Accelerator to execute
   * @param stream The file stream containing kernel source code
   */
  Program(std::shared_ptr<Accelerator> acc, std::istream &stream)
      : accelerator(std::move(acc)),
        src(std::istreambuf_iterator<char>(stream), {}) {
  }

  std::shared_ptr<IR> getIR() { return xaccIR; }

  /**
   * Execute the compilation mechanism on the provided program
   * source kernel code to produce XACC IR that can be executed
   * on the attached Accelerator.
   */
  virtual void build();
  virtual void build(const std::string &compiler);

  /**
   * Return the number of Kernels in this Program's XACC IR.
   *
   * @return nKernels The number of kernels.
   */
  const int nKernels() {
    if (!xaccIR)
      build();

    return xaccIR->getKernels().size();
  }

  /**
   * Return an executable version of the kernel
   * referenced by the kernelName string.
   *
   * @param name The name of the kernel
   * @return kernel The Kernel represented by kernelName
   */
  template <typename... RuntimeArgs>
  auto getKernel(const std::string &kernelName) -> Kernel<RuntimeArgs...> {

    // Build the kernel with the appropriate compiler
    if (!xaccIR) {
      build();
    }

    return Kernel<RuntimeArgs...>(accelerator, xaccIR->getKernel(kernelName));
  }

  /**
   * Return all Kernels that have sizeof...(RuntimeArgs)
   * InstructionParameters.
   * @return kernels Kernels with sizeof...(RuntimeArgs) Parameters.
   */
  template <typename... RuntimeArgs>
  auto getKernels() -> KernelList<RuntimeArgs...> {
    if (!xaccIR) {
      build();
    }

    KernelList<RuntimeArgs...> kernels(accelerator);
    for (auto k : xaccIR->getKernels()) {
      if (k->nParameters() == (sizeof...(RuntimeArgs))) {
        kernels.push_back(Kernel<RuntimeArgs...>(accelerator, k));
      }
    }

    return kernels;
  }

  /**
   * Return the Kernel range given by the provided indices.
   * @param beginIdx The beginning index for the range
   * @param endIdx The endindex for the range
   * @return kernels The kernels in the given range.
   */
  template <typename... RuntimeArgs>
  auto getKernels(const int beginIdx, const int endIdx)
      -> KernelList<RuntimeArgs...> {
    KernelList<RuntimeArgs...> kernels(accelerator);
    if (!xaccIR) {
      build();
    }

    for (int i = beginIdx; i < endIdx; i++) {
      if (xaccIR->getKernels()[i]->nParameters() == (sizeof...(RuntimeArgs))) {
        kernels.push_back(
            Kernel<RuntimeArgs...>(accelerator, xaccIR->getKernels()[i]));
      }
    }

    return kernels;
  }

  /**
   * Return a list of Kernels with now runtime type information
   * for the arguments. Any arguments must be provided
   * as InstructionParameters.
   */
  auto getRuntimeKernels() -> KernelList<> {
    if (!xaccIR) {
      build();
    }

    KernelList<> kernels(accelerator);
    for (auto k : xaccIR->getKernels()) {
      kernels.push_back(Kernel<>(accelerator, k));
    }

    return kernels;
  }

  virtual ~Program() {}
};

} // namespace xacc
#endif