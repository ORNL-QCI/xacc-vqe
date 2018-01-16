/***********************************************************************************
 * Copyright (c) 2016, UT-Battelle
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the xacc nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Contributors:
 *   Initial API and implementation - Alex McCaskey
 *
 **********************************************************************************/
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE FCIDumpPreprocessorTester

#include <boost/test/included/unit_test.hpp>
#include "FCIDumpPreprocessor.hpp"
#include "XACC.hpp"
#include <boost/mpi.hpp>

using namespace xacc::vqe;

using namespace boost;

// Create Global MPI Fixture to initialize MPI environment
struct MPIFixture {
	MPIFixture() :
			env(new mpi::environment()) {
		BOOST_TEST_MESSAGE("Setting up MPI Environment");
	}
	~MPIFixture() {
		BOOST_TEST_MESSAGE("Finalizing MPI Environment");
	}

	std::shared_ptr<mpi::environment> env;
};

// Make that fixture globals
BOOST_GLOBAL_FIXTURE(MPIFixture);

class FakeAcc: public xacc::Accelerator {
public:

	virtual std::shared_ptr<xacc::AcceleratorGraph> getAcceleratorConnectivity() {
		return std::make_shared<xacc::AcceleratorGraph>();
	}

	virtual xacc::AcceleratorType getType() {
		return xacc::AcceleratorType::qpu_gate;
	}

	virtual bool isValidBufferSize(const int nBits) {
		return true;
	}

	virtual std::vector<std::shared_ptr<xacc::IRTransformation>> getIRTransformations() {

	}
	;

	virtual std::shared_ptr<xacc::AcceleratorBuffer> createBuffer(
			const std::string& varId) {
		auto b = std::make_shared<xacc::AcceleratorBuffer>(varId, 1);
		storeBuffer(varId, b);
		return b;
	}

	virtual std::vector<std::string> getAllocatedBufferNames() {
		std::vector<std::string> names;
		names.push_back("hello");
		return names;
	}

	virtual std::shared_ptr<xacc::AcceleratorBuffer> getBuffer(
			const std::string& varid) {
		return std::make_shared<xacc::AcceleratorBuffer>("hello", 1);
	}

	virtual void initialize() {

	}
	/**
	 * Execute the provided XACC IR Function on the provided AcceleratorBuffer.
	 *
	 * @param buffer The buffer of bits this Accelerator should operate on.
	 * @param function The kernel to execute.
	 */
	virtual void execute(std::shared_ptr<xacc::AcceleratorBuffer> buffer,
			const std::shared_ptr<xacc::Function> function) {
	}

	virtual std::vector<std::shared_ptr<xacc::AcceleratorBuffer>> execute(
			std::shared_ptr<xacc::AcceleratorBuffer> buffer,
			const std::vector<std::shared_ptr<xacc::Function>> function) {
	}

	virtual const std::string name() const {
		return "fake";
	}

	virtual const std::string description() const {
		return "....";
	}

	/**
	 * Create, store, and return an AcceleratorBuffer with the given
	 * variable id string and of the given number of bits.
	 * The string id serves as a unique identifier
	 * for future lookups and reuse of the AcceleratorBuffer.
	 *
	 * @param varId The variable name of the created buffer
	 * @param size The number of bits in the created buffer
	 * @return buffer The buffer instance created.
	 */
	virtual std::shared_ptr<xacc::AcceleratorBuffer> createBuffer(
			const std::string& varId, const int size) {
		auto b = std::make_shared<xacc::AcceleratorBuffer>(varId, size);
		storeBuffer(varId, b);
		return b;
	}
};

class FakeCompiler: public xacc::Compiler {
public:
	virtual std::shared_ptr<xacc::IR> compile(const std::string& src,
			std::shared_ptr<xacc::Accelerator> acc) {

	}

	/**
	 * This method is to be implemented by derived Compilers
	 * and is in charge of executing the compilation mechanism
	 * on the provided source string.
	 * @param src
	 * @return
	 */
	virtual std::shared_ptr<xacc::IR> compile(const std::string& src) {

	}

	/**
	 * This method is to be implemented by derived Compilers and
	 * is in charge of taking the provided Function IR and converting
	 * it to source code in this Compiler's language.
	 *
	 * @param function The XACC IR Function to translate
	 * @return src The source code as a string
	 */
	virtual const std::string translate(const std::string& bufferVariable,
			std::shared_ptr<xacc::Function> function) {

	}

	/**
	 * Return the name of this Compiler
	 * @return name Compiler name
	 */
	virtual const std::string getName() {
		return "fake";
	}

	virtual const std::string name() const {
		return "fake";
	}

	/**
	 * Return the description of this instance
	 * @return description The description of this object.
	 */
	virtual const std::string description() const {
		return "";
	}

};

BOOST_AUTO_TEST_CASE(checkH2Process) {

	const std::string h2FciDump =
			R"h2FciDump(&FCI NORB=  2,NELEC=  2,MS2= 0,
  ORBSYM=1,5,
  ISYM=0,
  &END
  0.6744931033260081E+00   1   1   1   1  
  0.6634720448605567E+00   2   2   1   1  
  0.6973979494693358E+00   2   2   2   2  
  0.1812875358123322E+00   2   1   2   1  
 -0.1252477303982147E+01   1   1   0   0  
 -0.4759344611440753E+00   2   2   0   0  
  0.7137758743754461E+00   0   0   0   0)h2FciDump";

	auto compiler = std::make_shared<FakeCompiler>();
	auto acc = std::make_shared<FakeAcc>();

	FCIDumpPreprocessor processor;

	auto newKernelCode = processor.process(h2FciDump, compiler, acc);

	std::cout << "NEWKERNEL:\n" << newKernelCode << "\n";

	std::string expected =
			R"expected(__qpu__ kernel() {
   0.7137758743754461
   -1.252477303982147 0 1 0 0
   0.337246551663004 0 1 1 1 1 0 0 0
   0.0906437679061661 0 1 1 1 3 0 2 0
   0.0906437679061661 0 1 2 1 0 0 2 0
   0.3317360224302783 0 1 2 1 2 0 0 0
   0.0906437679061661 0 1 3 1 1 0 2 0
   0.3317360224302783 0 1 3 1 3 0 0 0
   0.337246551663004 1 1 0 1 0 0 1 0
   0.0906437679061661 1 1 0 1 2 0 3 0
   -1.252477303982147 1 1 1 0
   0.0906437679061661 1 1 2 1 0 0 3 0
   0.3317360224302783 1 1 2 1 2 0 1 0
   0.0906437679061661 1 1 3 1 1 0 3 0
   0.3317360224302783 1 1 3 1 3 0 1 0
   0.3317360224302783 2 1 0 1 0 0 2 0
   0.0906437679061661 2 1 0 1 2 0 0 0
   0.3317360224302783 2 1 1 1 1 0 2 0
   0.0906437679061661 2 1 1 1 3 0 0 0
   -0.4759344611440753 2 1 2 0
   0.0906437679061661 2 1 3 1 1 0 0 0
   0.3486989747346679 2 1 3 1 3 0 2 0
   0.3317360224302783 3 1 0 1 0 0 3 0
   0.0906437679061661 3 1 0 1 2 0 1 0
   0.3317360224302783 3 1 1 1 1 0 3 0
   0.0906437679061661 3 1 1 1 3 0 1 0
   0.0906437679061661 3 1 2 1 0 0 1 0
   0.3486989747346679 3 1 2 1 2 0 3 0
   -0.4759344611440753 3 1 3 0
}
)expected";

	BOOST_VERIFY(expected == newKernelCode);
}
