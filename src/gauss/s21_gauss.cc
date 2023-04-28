#include "sfleta_gauss.h"

namespace sfleta {
void Gauss::LoadMatrixFromFile(std::string filename) {
  std::fstream fs;
  fs.open(filename, std::fstream::in | std::fstream::out);
  if (!fs.is_open())
    throw std::invalid_argument(
        "Ошибка при открытии файла c месторасположением: " + filename);

  fs >> size_;
  matrix_.resize(size_, Vector(size_ + 1));

  for (auto &rows : matrix_) {
    fs >> rows;
    if (fs.fail())
      throw std::invalid_argument("Fail read data for the matrix!");
  }

  fs.close();

  count_thread_ = std::min(std::thread::hardware_concurrency(),
                           static_cast<unsigned int>(std::log2(size_)));
}

void Gauss::GaussianAlgorithm(launch launch) {
  gauss_mat_ = matrix_;
  result_.clear();

  for (int n = 0; n + shift_ < size_;) {
    if (IsNullColumn(n)) continue;

    gauss_mat_[n] /= gauss_mat_[n][n + shift_];
    if (launch == launch::ASYNC) {
      AsyncTriangulation(n);
    } else {
      Triangulation(n, 0, size_);
    }
    ++n;
  }
  ReverseStroke();
}

void Gauss::ReverseStroke() {
  for (int i = 0; i < size_; ++i) {
    double res = 0;
    for (int j = 0; j < size_; ++j) {
      res += matrix_[i][j] * gauss_mat_[j][size_];
    }
    if (abs(res - matrix_[i][size_]) > kEps) return;
  }
  CreateResultVector();
}

void Gauss::CreateResultVector() {
  std::transform(gauss_mat_.begin(), gauss_mat_.end(),
                 std::back_inserter(result_),
                 [this](auto rows) { return rows[size_]; });
}

void Gauss::AsyncTriangulation(int n) {
  std::vector<std::thread> vector_thread;
  for (int i = 0; i < count_thread_; ++i) {
    float div = static_cast<float>(size_) / count_thread_;
    vector_thread.push_back(
        std::thread(&Gauss::Triangulation, this, n, i * div, (i + 1) * div));
  }
  for (auto &thread : vector_thread) {
    thread.join();
  }
}

void Gauss::Triangulation(int n, int start, int end) {
  for (int i = start; i < end; ++i) {
    if (i == n) continue;
    double factor = gauss_mat_[i][n + shift_];
    gauss_mat_[i] -= std::move(gauss_mat_[n] * factor);
  }
}

bool Gauss::SwapRows(int n) {
  for (auto i = n; i < size_; ++i) {
    if (gauss_mat_[i][n + shift_]) {
      gauss_mat_[i].swap(gauss_mat_[n]);
      return true;
    }
  }
  ++shift_;
  return false;
}

bool Gauss::IsNullColumn(int n) {
  return !gauss_mat_[n][n + shift_] && !SwapRows(n);
}

void Gauss::PrintResult() {
  if (result_.empty()) {
    std::cout << "СЛАУ не имеет решений: " << std::endl;
  } else {
    std::cout << "Корни СЛАУ: " << std::endl;
    int index = 1;
    for (const auto &elem : result_) {
      std::cout << 'X' << index++ << " = ";
      std::cout << elem << '\t';
      if (!((index - 1) % 3)) std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
}  // namespace sfleta