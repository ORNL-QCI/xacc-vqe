#include "Program.hpp"
#include "XACC.hpp"
#include "Accelerator.hpp"
#include "IRTransformation.hpp"

#include <fstream>
#include "xacc_service.hpp"

namespace xacc {

std::vector<IRTransformation>
getAcceleratorIndependentTransformations(Accelerator::AcceleratorType &accType) {
  std::vector<IRTransformation> transformations;
  return transformations;
}

void Program::build(const std::string &compilerName) {
  // Get reference to the runtime options
  auto runtimeOptions = RuntimeOptions::instance();

  // Get the compiler that has been requested.
  auto compilerToRun =
      compilerName.empty() ? (*runtimeOptions)["compiler"] : compilerName;

  // Create the appropriate compiler
  compiler = xacc::getService<Compiler>(compilerToRun);

  // Make sure we got a valid
  if (!compiler) {
    xacc::error("Invalid Compiler requested in Program.build().\n");
  }

  XACCLogger::instance()->info("Executing " + compiler->name() + " compiler.");

  // Execute the compilation
  xaccIR = compiler->compile(src, accelerator);
  XACCLogger::instance()->info("Done executing " + compiler->name() +
                               " compiler.");

  // Validate the compilation
  if (!xaccIR) {
    XACCLogger::instance()->error("Bad source string or something.\n");
  }

  // Execute hardware dependent IR Transformations
  auto accTransforms = accelerator->getIRTransformations();
  for (auto t : accTransforms) {
    xaccIR = t->transform(xaccIR);
  }

  // Write the IR to file if the user requests it
  if (runtimeOptions->exists("persist-ir")) {
    auto fileStr = (*runtimeOptions)["persist-ir"];
    std::ofstream ostr(fileStr);
    xaccIR->persist(ostr);
  }
}

void Program::build() {
  build("");
  return;
}

} // namespace xacc