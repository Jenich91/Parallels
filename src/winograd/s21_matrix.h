#ifndef sfleta_MATRIX_H
#define sfleta_MATRIX_H
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

namespace sfleta {
class Matrix {
 public:
  using matrix_t = std::vector<std::vector<double> >;

  explicit Matrix(size_t rows, size_t cols);
  explicit Matrix(const std::string& filename);
  ~Matrix() {}

  void RandomFill();

  matrix_t getBuffer() const { return buffer_; }
  size_t getCols() const { return cols_; }
  size_t getRows() const { return rows_; }

  void mul_matrix(const Matrix& other);
  void Print();
  double& operator()(size_t, size_t);
  double operator()(size_t, size_t) const;

  Matrix& operator=(const Matrix& other);

 private:
  matrix_t buffer_;
  size_t rows_;
  size_t cols_;

  void Clear();
  void LoadFromFile(const std::string& filename);
};
}  // namespace sfleta
#endif  // sfleta_MATRIX_H
