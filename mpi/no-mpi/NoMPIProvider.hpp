#ifndef MPI_NO_MPI_NOMPIPROVIDER_HPP_
#define MPI_NO_MPI_NOMPIPROVIDER_HPP_

#include "MPIProvider.hpp"

namespace xacc {
namespace vqe {

class NullCommunicator : public Communicator {

public:

	NullCommunicator() {}

	virtual const int rank() {
		return 0;
	}

	virtual const int size() {
		return 1;
	}

	virtual void broadcast(std::vector<double>& d, const int r) {
		return;
	}

	virtual void broadcast(std::string& s, const int r) {
		return;
	}

	virtual void sumDoubles(double& myVal, double& result) {
		result = myVal;
	}

	virtual void sumInts(int& myVal, int& result) {
		result = myVal;
	}

	virtual ~NullCommunicator() {}

};

class NoMPIProvider : public MPIProvider {

protected:

	std::shared_ptr<Communicator> comm;
public:
	virtual void initialize() {
		comm = std::make_shared<NullCommunicator>();
	}

	virtual void initialize(int argc, char** argv) {
		comm = std::make_shared<NullCommunicator>();
	}

	virtual std::shared_ptr<Communicator> getCommunicator() {
		return comm;
	}

	virtual const std::string name() const {
		return "no-mpi";
	}

	virtual const std::string description() const {
		return "";
	}

	virtual ~NoMPIProvider() {}
};

}
}
#endif
