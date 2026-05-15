#pragma once

#include <vector>
#include <utility>

#include "task/include/task.hpp"
#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"

namespace kutergin_a_multidim_trapezoid { 

class KuterginAMultidimTrapezoidSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }

  explicit KuterginAMultidimTrapezoidSTL(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  double CalculateChunkSum(size_t start_idx, size_t end_idx, const std::vector<double>& h);

  double res_{0.0};
};

}  // namespace kutergin_a_multidim_trapezoid