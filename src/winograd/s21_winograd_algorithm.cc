#include "sfleta_winograd_algorithm.h"

namespace sfleta {
Matrix *Winograd::winograd_algo(const Matrix &matrix_, const Matrix &other) {
  if (matrix_.getRows() != other.getCols() ||
      matrix_.getCols() != other.getRows()) {
    throw std::invalid_argument(
        "The rows/columns count first matrix is not equal rows/columns of the "
        "second matrix");
  }

  size_t a = matrix_.getRows();
  size_t b = matrix_.getCols();
  size_t c = other.getCols();
  size_t d = matrix_.getCols() / 2;

  std::vector<double> rowsResult(a);
  std::vector<double> columnsResult(c);

  const auto &A = matrix_;
  const auto &B = other;
  Matrix &R = *new Matrix(a, c);

  // алгоритм корректно работает с матрицами где b>1
  if (b == 1) {
    R = matrix_;
    R.mul_matrix(other);
    return &R;
  }

  for (size_t i = 0; i < a; ++i) {
    rowsResult[i] = A(i, 0) * A(i, 1);
    for (size_t j = 1; j < d; ++j) {
      rowsResult[i] += A(i, 2 * j) * A(i, 2 * j + 1);
    }
  }

  for (size_t i = 0; i < c; ++i) {
    columnsResult[i] = B(0, i) * B(1, i);
    for (size_t j = 1; j < d; ++j) {
      columnsResult[i] += B(2 * j, i) * B(2 * j + 1, i);
    }
  }

  for (size_t i = 0; i < a; ++i) {
    for (size_t j = 0; j < c; ++j) {
      R(i, j) = -rowsResult[i] - columnsResult[j];
      for (size_t k = 0; k < d; ++k) {
        R(i, j) +=
            (A(i, 2 * k) + B(2 * k + 1, j)) * (A(i, 2 * k + 1) + B(2 * k, j));
      }
    }
  }

  //  прибавление членов в случае нечетной общей размерности
  if (b % 2 != 0) {
    for (size_t i = 0; i < a; ++i) {
      for (size_t j = 0; j < c; ++j) {
        R(i, j) += A(i, b - 1) * B(b - 1, j);
      }
    }
  }

  return &R;
}

Matrix *Winograd::winograd_algo_multithread(const Matrix &matrix_,
                                            const Matrix &other,
                                            size_t thread_amount) {
  if (matrix_.getRows() != other.getCols() ||
      matrix_.getCols() != other.getRows()) {
    throw std::invalid_argument(
        "The rows/columns count first matrix is not equal rows/columns of the "
        "second matrix");
  }

  if (thread_amount == 0) {
    throw std::invalid_argument("Invalid threads amount");
  }
  if (thread_amount == 1) {
    return winograd_algo(matrix_, other);
  }

  size_t a = matrix_.getRows();
  size_t b = matrix_.getCols();
  size_t c = other.getCols();
  size_t d = matrix_.getCols() / 2;

  auto &A = matrix_;
  auto &B = other;
  Matrix &R = *new Matrix(a, c);

  // алгоритм корректно работает с матрицами где b>1
  if (b == 1) {
    R = matrix_;
    R.mul_matrix(other);
    return &R;
  }

  std::vector<double> rowsResult(a);
  std::vector<double> columnsResult(c);

  std::thread thread1([&]() { PreCalcRowsResult(A, rowsResult, 0, a, d); });
  std::thread thread2(
      [&]() { PreCalcColumnsResult(B, columnsResult, 0, c, d); });
  thread1.join();
  thread2.join();

  std::vector<std::thread> threads;
  size_t count = 0;
  for (size_t from = 0, step = a / thread_amount, to = step;
       count < thread_amount; from += step, to += step) {
    if (++count == thread_amount && to < a) {
      to = a;
    }

    threads.push_back(std::thread(&Winograd::PreCalcFinal, std::ref(rowsResult),
                                  std::ref(columnsResult), std::ref(A),
                                  std::ref(B), std::ref(R), from, to, b, c, d));
  }
  for (size_t k = 0; k < count; ++k) {
    threads[k].join();
  }

  return &R;
}

void Winograd::PreCalcFinal(const std::vector<double> &rowsResult,
                            const std::vector<double> &columnsResult,
                            const Matrix &A, const Matrix &B, Matrix &R,
                            size_t i, size_t a, size_t b, size_t c, size_t d) {
  for (; i < a; ++i) {
    for (size_t j = 0; j < c; ++j) {
      R(i, j) = -rowsResult[i] - columnsResult[j];
      if (b % 2 != 0) {
        R(i, j) += A(i, b - 1) * B(b - 1, j);
      }
      for (size_t k = 0; k < d; ++k) {
        R(i, j) +=
            (A(i, 2 * k) + B(2 * k + 1, j)) * (A(i, 2 * k + 1) + B(2 * k, j));
      }
    }
  }
}

void Winograd::PreCalcRowsResult(const Matrix &A,
                                 std::vector<double> &rowsResult, size_t i,
                                 size_t a, size_t d) {
  for (; i < a; ++i) {
    rowsResult[i] = A(i, 0) * A(i, 1);
    for (size_t j = 1; j < d; ++j) {
      rowsResult[i] += A(i, 2 * j) * A(i, 2 * j + 1);
    }
  }
}

void Winograd::PreCalcColumnsResult(const Matrix &B,
                                    std::vector<double> &columnsResult,
                                    size_t i, size_t c, size_t d) {
  for (; i < c; ++i) {
    columnsResult[i] = B(0, i) * B(1, i);
    for (size_t j = 1; j < d; ++j) {
      columnsResult[i] += B(2 * j, i) * B(2 * j + 1, i);
    }
  }
}
}  // namespace sfleta