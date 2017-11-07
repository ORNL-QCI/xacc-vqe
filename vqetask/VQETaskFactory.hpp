/*
 * VQETask.hpp
 *
 *  Created on: Nov 7, 2017
 *      Author: aqw
 */

#ifndef VQETASK_VQETASK_HPP_
#define VQETASK_VQETASK_HPP_

#include "ComputeEnergyVQETask.hpp"

namespace xacc {
namespace vqe {

class VQETaskFactory {

public:

	virtual std::shared_ptr<VQETask> createTask(const std::string& taskType, std::shared_ptr<VQEProgram> program) {
		if (taskType == "ComputeEnergy") {
			return std::make_shared<ComputeEnergyVQETask>(program);
		}
	}

};

}
}

#endif /* IR_VQETASK_HPP_ */
