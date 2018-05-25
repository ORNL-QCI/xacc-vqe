#ifndef IR_VQETASK_HPP_
#define IR_VQETASK_HPP_

#include "OptionsProvider.hpp"
#include <Eigen/Dense>
#include "VQEProgram.hpp"

namespace xacc {
namespace vqe {

class VQETaskResult {

protected:

	std::string _fileName;

	//std::ofstream _output;

public:
	VQETaskResult() {}

	VQETaskResult(const std::string& fileName) : 
		_fileName(fileName) {}

	void persist() {
		if(!_fileName.empty()) { // && _output.is_open()) {
			std::stringstream ss;
			// check if its been written to already
			if (std::ifstream(_fileName).peek() == std::ifstream::traits_type::eof()) {
				// write header
				for (int i = 0; i < angles.size(); i++) ss << "t" << i << ",";
				for (auto& kv : readoutErrorProbabilities) ss << kv.first << ",";
				for (auto& kv : expVals) ss << kv.first << ",";
				ss << "E\n";
			}

			// write angles
			for (int i = 0; i < angles.size(); i++) ss << angles(i) << ",";
			// write probabilities
			for (auto& kv : readoutErrorProbabilities) ss << kv.second << ",";
			// write exp vals
			for (auto& kv : expVals) ss << kv.second << ",";
			// write energy
			ss << energy << "\n";

			std::ofstream _output(_fileName, std::ofstream::app);
			_output << ss.str();
			_output.close();
		}
	}

	std::map<std::string, double> expVals;

	double energy = 0.0;

	Eigen::VectorXd angles;

	int nQpuCalls = 0;

	int vqeIterations = 0;

	std::map<std::string, double> readoutErrorProbabilities;
};

class VQETask : public xacc::Identifiable, public OptionsProvider {

public:

	VQETask() {}

	VQETask(std::shared_ptr<VQEProgram> prog) : program(prog) {
	}

	virtual VQETaskResult execute(Eigen::VectorXd parameters) = 0;

	virtual void setVQEProgram(std::shared_ptr<VQEProgram> p) {
		program = p;
	}

	/**
	 * Return an empty options_description, this is for
	 * subclasses to implement.
	 */
	virtual std::shared_ptr<options_description> getOptions() {
		return std::make_shared<options_description>();
	}

	/**
	 * Given user-input command line options, perform
	 * some operation. Returns true if runtime should exit,
	 * false otherwise.
	 *
	 * @param map The mapping of options to values
	 * @return exit True if exit, false otherwise
	 */
	virtual bool handleOptions(variables_map& map) {
		return false;
	}

	virtual ~VQETask() {}

protected:

	std::shared_ptr<VQEProgram> program;
};

}
}

#endif
