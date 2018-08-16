#include "UCCSD.hpp"
#include "GateFunction.hpp"
#include "FermionToSpinTransformation.hpp"
#include "CommutingSetGenerator.hpp"
#include <boost/math/constants/constants.hpp>

using namespace xacc::quantum;

namespace xacc {
namespace vqe {

std::shared_ptr<Function> UCCSD::generate(
			std::map<std::string, InstructionParameter> parameters) {

    if (!parameters.count("n-electrons") && !parameters.count("n_electrons")) {
        xacc::error("Invalid mapping of parameters for UCCSD generator, missing n-electrons key.");
    }

    if (!parameters.count("n-qubits") && !parameters.count("n_qubits")) {
        xacc::error("Invalid mapping of parameters for UCCSD generator, missing n-qubits key.");
    }

    bool hasUnderscore = false;

    if (parameters.count("n_electrons")) {
        hasUnderscore = true;
        if (!parameters.count("n_qubits")) xacc::error("UCCSD Generator missing n_qubits key.");
    }

    if (parameters.count("n_qubits")) {
        hasUnderscore = true;
        if (!parameters.count("n_electrons")) xacc::error("UCCSD Generator missing n_electrons key.");
    }

    std::vector<InstructionParameter> params;
    if (hasUnderscore) {
        xacc::setOption("n-electrons",boost::lexical_cast<std::string>(parameters["n_electrons"]));
        params.push_back(parameters["n_electrons"]);
        params.push_back(parameters["n_qubits"]);
    } else {
        params.push_back(parameters["n-electrons"]);
        params.push_back(parameters["n-qubits"]);
    }

    if (parameters.size() > 2) {
        for (auto& kv : parameters) {
            if (!boost::contains(kv.first, "electrons") && !boost::contains(kv.first, "qubits")) {
                params.push_back(kv.second);
            }
        }
    }
    
    return generate(nullptr, params);
}

std::shared_ptr<Function> UCCSD::generate(
		std::shared_ptr<AcceleratorBuffer> buffer,
		std::vector<InstructionParameter> parameters) {

    xacc::info("Running UCCSD Generator.");
    std::vector<xacc::InstructionParameter> variables;

    int nQubits = 0;
    int nElectrons = 0;
    if (parameters.empty()) {
    	if (!xacc::optionExists("n-electrons")) {
	    	xacc::error("To use this UCCSD State Prep IRGenerator, you "
		    		"must specify the number of electrons.");
	    }

	    if (!xacc::optionExists("n-qubits")) {
		    xacc::error("To use this UCCSD State Prep IRGenerator, you "
			    	"must specify the number of qubits.");
	    }

	    nQubits = std::stoi(xacc::getOption("n-qubits"));
	    nElectrons = std::stoi(xacc::getOption("n-electrons"));

    } else {
        if (parameters.size() < 2) xacc::error("Invalid input parameters for UCCSD generator.");

        nElectrons = boost::get<int>(parameters[0]);
        nQubits = boost::get<int>(parameters[1]);

        if (parameters.size() > 2) {
            for (int i = 2; i < parameters.size(); i++) {
                variables.push_back(parameters[i]);
            }
        }
    } 
    
    xacc::info("UCCSD Generator (nqubits,nelectrons) = " + std::to_string(nQubits)+", " + std::to_string(nElectrons) +".");

	// Compute the number of parameters
	auto _nOccupied = (int) std::ceil(nElectrons / 2.0);
	auto _nVirtual = nQubits / 2 - _nOccupied;
	auto nSingle = _nOccupied * _nVirtual;
	auto nDouble = nSingle * (nSingle+1) / 2; 
	auto _nParameters = nSingle + nDouble;

	std::vector<std::string> params;

    if (variables.empty()) {
	    for (int i = 0; i < _nParameters; i++) {
            params.push_back("theta" + std::to_string(i));
		    variables.push_back(InstructionParameter("theta" + std::to_string(i)));
        }
	} else {
        for (int i = 0; i < _nParameters; i++) {
		    params.push_back(boost::get<std::string>(variables[i]));
	    }
    }

    auto slice = [](const std::vector<std::string>& v, int start=0, int end=-1) {
        int oldlen = v.size();
        int newlen;
        if (end == -1 or end >= oldlen){
            newlen = oldlen-start;
        } else {
            newlen = end-start;
        }
        std::vector<std::string> nv(newlen);
        for (int i=0; i<newlen; i++) {
            nv[i] = v[start+i];
        }
        return nv;
    };
    auto singleParams = slice(params, 0, nSingle);
    auto doubleParams1 = slice(params, nSingle, 2*nSingle);
    auto doubleParams2 = slice(params, 2*nSingle);
    std::vector<std::function<int(int)>> fs{[](int i) {return 2*i;}, [](int i) {return 2*i+1;}};
    
    using OpType = std::vector<std::pair<int,int>>;
	auto kernel = std::make_shared<FermionKernel>("fermiUCCSD");
    int count = 0;
    for (int i = 0; i < _nVirtual; i++) {
        for (int j = 0; j < _nOccupied; j++) {
            auto vs = _nOccupied+i;
            auto os = j;
            for (int s = 0; s < 2; s++) {
                auto ti = fs[s];
                auto oi = fs[1-s];
                auto vt = ti(vs);
                auto vo = oi(vs);
                auto ot = ti(os);
                auto oo = oi(os);
                
                OpType op1{{vt,1},{ot,0}}, op2{{ot,1},{vt,0}};
                auto i1 = std::make_shared<FermionInstruction>(op1, singleParams[count]);
                auto i2 = std::make_shared<FermionInstruction>(op2, singleParams[count], std::complex<double>(-1.,0.));
                kernel->addInstruction(i1);
                kernel->addInstruction(i2);
                
                OpType op3{{vt,1},{ot,0}, {vo,1},{oo,0}}, op4{{oo,1}, {vo,0}, {ot,1},{vt,0}};
                auto i3 = std::make_shared<FermionInstruction>(op3, doubleParams1[count]);
                auto i4 = std::make_shared<FermionInstruction>(op4, doubleParams1[count], std::complex<double>(-1.,0.));

                kernel->addInstruction(i3);
                kernel->addInstruction(i4);

            }
            count++;
        }
    }

    count = 0;
    // routine for converting amplitudes for use by UCCSD
	std::vector< std::tuple<int, int> > tupleVec;
	for (int i = 0; i < _nVirtual; i++){
		for (int j = 0; j < _nOccupied; j++){
			tupleVec.push_back(std::make_tuple(i, j));
		}
	}
	// Combination lambda used to determine indices
	auto Combination = [=](std::vector< std::tuple<int, int> > t){
	std::vector< std::tuple<int, int, int, int> > comboVec;
		for (int i = 0; i < t.size(); i++){
			for (int j = i+1; j < t.size(); j++){
				const std::tuple<int, int, int, int> newTuple = std::tuple_cat(t[i], t[j]);
				comboVec.push_back(newTuple);
			}
		}
		return comboVec;
	};

	auto combineVec = Combination(tupleVec);
    for (auto i : combineVec) {
        auto p = std::get<0>(i);
        auto q = std::get<1>(i);
        auto r = std::get<2>(i);
        auto s = std::get<3>(i);

        auto vs1 = _nOccupied + p;
        auto os1 = q;
        auto vs2 = _nOccupied + r;
        auto os2 = s;

        for (int sa = 0; sa < 2; sa++) {
            for (int sb = 0; sb < 2; sb++) {
                auto ia = fs[sa];
                auto ib = fs[sb];

                auto v1a = ia(vs1);
                auto o1a = ia(os1);
                auto v2b = ib(vs2);
                auto o2b = ib(os2);
                
                OpType op5{{v1a,1},{o1a,0},{v2b,1},{o2b,0}}, op6{{o2b,1},{o1a,0},{v2b,1},{o2b,0}};
                auto i5 = std::make_shared<FermionInstruction>(op5, doubleParams2[count]);
                auto i6 = std::make_shared<FermionInstruction>(op6, doubleParams2[count], std::complex<double>(-1.,0.));

                kernel->addInstruction(i5);
                kernel->addInstruction(i6);
            }
        }
        count++;
    }
    
	// std::cout << "KERNEL: \n" << kernel->toString("") << "\n";
	// Create the FermionIR to pass to our transformation.
	auto fermionir = std::make_shared<FermionIR>();
	fermionir->addKernel(kernel);

	xacc::info("Done constructing UCCSD Fermion Operator.");
	xacc::info("Mapping UCCSD Fermion Operator to Spin. ");

	std::shared_ptr<FermionToSpinTransformation> transform;
	if (xacc::optionExists("fermion-transformation")) {
		auto transformStr = xacc::getOption("fermion-transformation");
		transform = xacc::getService<FermionToSpinTransformation>(
				transformStr);
	} else {
		transform = xacc::getService<FermionToSpinTransformation>(
				"jw");
	}

	auto compositeResult = transform->transform(*kernel.get());
//	auto resultsStr = compositeResult.toString();
//	boost::replace_all(resultsStr, "+", "+\n");

	// Create the Spin Hamiltonian
	auto transformedIR = compositeResult.toXACCIR();
	xacc::info("Done mapping UCCSD Fermion Operator to Spin.");

	std::unordered_map<std::string, Term> terms = compositeResult.getTerms();

	CommutingSetGenerator gen;
	auto commutingSets = gen.getCommutingSet(compositeResult, nQubits);
	auto pi = boost::math::constants::pi<double>();
	auto gateRegistry = xacc::getService<IRProvider>("gate");

	auto uccsdGateFunction = gateRegistry->createFunction("uccsdPrep",{},
			variables);


	// Perform Trotterization...
	for (auto s : commutingSets) {

		for (auto inst : s) {
			Term spinInst = inst;

			// Get the individual pauli terms
			auto termsMap = std::get<2>(spinInst);

			std::vector<std::pair<int,std::string>> terms;
			for (auto& kv : termsMap) {
				if (kv.second != "I" && !kv.second.empty()) {
					terms.push_back({kv.first, kv.second});
				}
			}
			// The largest qubit index is on the last term
			int largestQbitIdx = terms[terms.size() - 1].first;
			auto tempFunction = gateRegistry->createFunction("temp", {}, {});

			for (int i = 0; i < terms.size(); i++) {

				auto qbitIdx = terms[i].first;
				auto gateName = terms[i].second;

				if (i < terms.size() - 1) {
					auto cnot =
							gateRegistry->createInstruction(
									"CNOT",
									std::vector<int> { qbitIdx,
											terms[i + 1].first });
					tempFunction->addInstruction(cnot);
				}

				if (gateName == "X") {
					auto hadamard =
							gateRegistry->createInstruction(
									"H", std::vector<int> { qbitIdx });
					tempFunction->insertInstruction(0, hadamard);
				} else if (gateName == "Y") {
					auto rx =
							gateRegistry->createInstruction(
									"Rx", std::vector<int> { qbitIdx });
					InstructionParameter p(pi / 2.0);
					rx->setParameter(0, p);
					tempFunction->insertInstruction(0, rx);
				}

				// Add the Rotation for the last term
				if (i == terms.size() - 1) {
					// FIXME DONT FORGET DIVIDE BY 2
					std::stringstream ss;
					ss << 2*std::imag(std::get<0>(spinInst)) << " * "
							<< std::get<1>(spinInst);
					auto rz =
							gateRegistry->createInstruction(
									"Rz", std::vector<int> { qbitIdx });

					InstructionParameter p(ss.str());
					rz->setParameter(0, p);
					tempFunction->addInstruction(rz);

				}
			}

			int counter = tempFunction->nInstructions();
			// Add the instruction on the backend of the circuit
			for (int i = terms.size() - 1; i >= 0; i--) {

				auto qbitIdx = terms[i].first;
				auto gateName = terms[i].second;

				if (i < terms.size() - 1) {
					auto cnot =
							gateRegistry->createInstruction(
									"CNOT",
									std::vector<int> { qbitIdx,
											terms[i + 1].first });
					tempFunction->insertInstruction(counter, cnot);
					counter++;
				}

				if (gateName == "X") {
					auto hadamard =
							gateRegistry->createInstruction(
									"H", std::vector<int> { qbitIdx });
					tempFunction->addInstruction(hadamard);
				} else if (gateName == "Y") {
					auto rx =
							gateRegistry->createInstruction(
									"Rx", std::vector<int> { qbitIdx });
					InstructionParameter p(4 * pi - (pi /2.0));
					rx->setParameter(0, p);
					tempFunction->addInstruction(rx);
				}

			}
			// Add to the total UCCSD State Prep function
			for (auto inst : tempFunction->getInstructions()) {
				uccsdGateFunction->addInstruction(inst);
			}
		}
	}

	for (int i = nElectrons-1; i >= 0; i--) {
		auto xGate = gateRegistry->createInstruction(
				"X", std::vector<int>{i});
		uccsdGateFunction->insertInstruction(0,xGate);
	}

	return uccsdGateFunction;
}

}
}

