#pragma once
#include <chrono>
#include <iostream>

namespace sfleta {
class Timer {
 public:
  Timer() : start{std::chrono::high_resolution_clock::now()} {}
  ~Timer() {
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                                       start)
                     .count()
              << " microseconds" << std::endl;
  }

 private:
  std::chrono::_V2::high_resolution_clock::time_point start;
};
}  // namespace sfleta
