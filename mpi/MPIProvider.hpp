#ifndef MPI_MPIPROVIDER_HPP_
#define MPI_MPIPROVIDER_HPP_

#include "Identifiable.hpp"

namespace xacc {
namespace vqe {
class Communicator {
public:

	virtual const int rank() {
		return 0;
	}
	virtual const int size() {
		return 1;
	}
	virtual void broadcast(std::vector<double>& d, const int r) = 0;
	virtual void broadcast(std::string& d, const int r) = 0;

	virtual void sumDoubles(double& myVal, double& result) = 0;
	virtual void sumInts(int& myVal, int& result) = 0;

	virtual ~Communicator() {}

};
class MPIProvider : public xacc::Identifiable {

public:

	virtual void initialize(int argc, char** argv) = 0;
	virtual void initialize() = 0;

	virtual std::shared_ptr<Communicator> getCommunicator() = 0;

	virtual ~MPIProvider() {}

};

}
}


#endif
