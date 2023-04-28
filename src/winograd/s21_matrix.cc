#include "sfleta_matrix.h"

namespace sfleta {
Matrix::Matrix(size_t rows, size_t cols) {
  if (rows < 1 || cols < 1) {
    throw std::invalid_argument("Wrong matrix size");
  }

  rows_ = rows;
  cols_ = cols;

  buffer_.resize(rows, std::vector<double>(cols));
}

Matrix::Matrix(const std::string& filename) { LoadFromFile(filename); }

void Matrix::RandomFill() {
  static std::mt19937 dice(time(nullptr));
  std::for_each(buffer_.begin(), buffer_.end(), [](auto&& iv) {
    for (auto&& i : iv) {
      i = [](double fMin, double fMax) -> double {
        double f = (double)dice() / RAND_MAX;
        return fMin + f * (fMax - fMin);
      }(0, 100);
    }
  });
}

void Matrix::LoadFromFile(const std::string& filename) {
  std::fstream fs;
  fs.open(filename, std::fstream::in | std::fstream::out);
  if (!fs.is_open())
    throw std::invalid_argument(
        "Ошибка при открытии файла c месторасположением: " + filename);

  Clear();

  fs >> rows_;
  fs >> cols_;
  buffer_.resize(rows_, std::vector<double>(cols_));

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      fs >> buffer_[i][j];
      if (fs.fail()) {
        throw std::invalid_argument("Fail read data for the matrix!");
      }
    }
  }
  fs.close();
}

void Matrix::Clear() {
  if (buffer_.size()) {
    for (auto& row : buffer_) row.clear();
    buffer_.clear();
  }
  rows_ = 0;
  cols_ = 0;
}

void Matrix::Print() {
  std::for_each(buffer_.begin(), buffer_.end(), [](const auto& iv) {
    for (const auto& i : iv) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  });
  std::cout << std::endl;
}

double& Matrix::operator()(size_t rows, size_t cols) {
  if (rows >= rows_ || cols >= cols_) {
    throw std::out_of_range("Array index out of bound");
  }
  return buffer_[rows][cols];
}

double Matrix::operator()(size_t rows, size_t cols) const {
  if (rows >= rows_ || cols >= cols_) {
    throw std::out_of_range("Array index out of bound");
  }
  return buffer_[rows][cols];
}

void Matrix::mul_matrix(const Matrix& other) {
  if (rows_ != other.cols_ || cols_ != other.rows_) {
    throw std::invalid_argument(
        "The rows/columns count first matrix is not equal rows/columns of the "
        "second matrix");
  }

  Matrix tempMatrix(rows_, other.cols_);
  for (size_t i = 0; i < rows_; i++) {
    for (size_t j = 0; j < other.cols_; j++) {
      for (size_t k = 0; k < cols_; k++) {
        tempMatrix.buffer_[i][j] += buffer_[i][k] * other.buffer_[k][j];
      }
    }
  }

  *this = tempMatrix;
}

Matrix& Matrix::operator=(const Matrix& other) {
  if (this == &other) {
    return *this;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  buffer_ = other.getBuffer();

  return *this;
}

}  // namespace sfleta
