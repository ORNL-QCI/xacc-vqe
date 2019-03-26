/*******************************************************************************
 * Copyright (c) 2018 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Alexander J. McCaskey - initial API and implementation
 *******************************************************************************/
#ifndef XACC_VQERESTARTDECORATOR_HPP_
#define XACC_VQERESTARTDECORATOR_HPP_

#include "AcceleratorDecorator.hpp"

namespace xacc {

namespace vqe {

class VQERestartDecorator : public AcceleratorDecorator {
public:
  void initialize() override;

  void execute(std::shared_ptr<AcceleratorBuffer> buffer,
               const std::shared_ptr<Function> function) override;

  std::vector<std::shared_ptr<AcceleratorBuffer>>
  execute(std::shared_ptr<AcceleratorBuffer> buffer,
          const std::vector<std::shared_ptr<Function>> functions) override;

  const std::string name() const override { return "vqe-restart"; }
  const std::string description() const override { return ""; }

  OptionPairs getOptions() override {
    OptionPairs desc {{"vqe-restart-file",
                        ""}};
    return desc;
  }
  ~VQERestartDecorator() override {}

private:

  std::shared_ptr<AcceleratorBuffer> loadedBuffer;
  bool initialized = false;
  std::queue<std::vector<std::shared_ptr<AcceleratorBuffer>>> queuedChildren;
};

} // namespace vqe
} // namespace xacc
#endif
