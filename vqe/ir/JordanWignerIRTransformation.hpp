#ifndef VQE_IR_JORDANWIGNERIRTRANSFORMATION_HPP_
#define VQE_IR_JORDANWIGNERIRTRANSFORMATION_HPP_

#include "IRTransformation.hpp"
#include "FermionKernel.hpp"
#include "FermionIR.hpp"

namespace xacc {

namespace vqe {

class JordanWignerIRTransformation: public xacc::IRTransformation {

public:

	virtual std::shared_ptr<IR> transform(std::shared_ptr<IR> ir);
};

}

}

#endif
