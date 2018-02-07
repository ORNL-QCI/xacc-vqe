#ifndef MPI_BOOST_MPI_BOOSTMPIPROVIDER_HPP_
#define MPI_BOOST_MPI_BOOSTMPIPROVIDER_HPP_


#include "MPIProvider.hpp"
#include <boost/mpi.hpp>

namespace xacc {
namespace vqe {

class BoostCommunicator : public Communicator {

protected:
	boost::mpi::communicator comm;
public:

	BoostCommunicator(boost::mpi::communicator& c) :comm(c) {}

	virtual const int rank() {
		return comm.rank();
	}

	virtual const int size() {
		return comm.size();
	}

	virtual void broadcast(std::vector<double>& d, const int r) {
		boost::mpi::broadcast(comm, d, r);
	}

	virtual void broadcast(std::string& s, const int r) {
		boost::mpi::broadcast(comm, s, r);
	}

	virtual void sumDoubles(double& myVal, double& result) {
		boost::mpi::all_reduce(comm, myVal, result, std::plus<double>());
	}

	virtual void sumInts(int& myVal, int& result) {
		boost::mpi::all_reduce(comm, myVal, result, std::plus<int>());
	}

	virtual ~BoostCommunicator() {}

};

class BoostMPIProvider : public MPIProvider {

protected:

	std::shared_ptr<boost::mpi::environment> env;

	std::shared_ptr<Communicator> comm;
public:

	virtual void initialize() {
		if (!env) env = std::make_shared<boost::mpi::environment>();
		boost::mpi::communicator c;
		comm = std::make_shared<BoostCommunicator>(c);
	}

	virtual void initialize(int argc, char** argv) {
		env = std::make_shared<boost::mpi::environment>(argc, argv);
		boost::mpi::communicator c;
		comm = std::make_shared<BoostCommunicator>(c);
	}

	virtual std::shared_ptr<Communicator> getCommunicator() {
		return comm;
	}

	virtual const std::string name() const {
		return "boost-mpi";
	}

	virtual const std::string description() const {
		return "";
	}

	virtual ~BoostMPIProvider() {}
};

}
}
#endif
