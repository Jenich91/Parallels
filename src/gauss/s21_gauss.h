#pragma once

#include <math.h>

#include <thread>

#include "sfleta_vector.h"

namespace sfleta {
static const constexpr double kEps = 1e-6;
class Gauss {
 public:
  enum launch { ASYNC = true, DEFFERED = false };

  Gauss() : size_(0), shift_(0), count_thread_(0){};
  ~Gauss() = default;
  Gauss(const Gauss&) = default;

  void GaussianAlgorithm(launch launch = DEFFERED);
  Vector& GetResult() { return result_; };
  void LoadMatrixFromFile(std::string);
  void PrintResult();

 private:
  std::vector<Vector> matrix_;
  std::vector<Vector> gauss_mat_;
  Vector result_;

  int size_;
  int shift_;
  int count_thread_;

  void Triangulation(int, int, int);
  void AsyncTriangulation(int);
  void CreateResultVector();
  bool IsNullColumn(int);
  void ReverseStroke();
  bool SwapRows(int);
};
}  // namespace sfleta
