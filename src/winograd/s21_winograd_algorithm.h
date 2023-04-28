#ifndef A3_PARALLELS_0_WINOGRAD_ALGORITHM_H
#define A3_PARALLELS_0_WINOGRAD_ALGORITHM_H
#include "sfleta_matrix.h"

namespace sfleta {
class Winograd {
 public:
  static Matrix* winograd_algo(const Matrix&, const Matrix&);
  static Matrix* winograd_algo_multithread(
      const Matrix&, const Matrix&,
      size_t thread_amount = std::thread::hardware_concurrency());

 private:
  Winograd() {}

  static void PreCalcRowsResult(const Matrix& A,
                                std::vector<double>& rowsResult, size_t i,
                                size_t a, size_t d);
  static void PreCalcColumnsResult(const Matrix& B,
                                   std::vector<double>& columnsResult, size_t i,
                                   size_t c, size_t d);
  static void PreCalcFinal(const std::vector<double>& rowsResult,
                           const std::vector<double>& columnsResult,
                           const Matrix& A, const Matrix& B, Matrix& R,
                           size_t i, size_t a, size_t b, size_t c, size_t d);
};
}  // namespace sfleta
#endif  // A3_PARALLELS_0_WINOGRAD_ALGORITHM_H
