#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace sfleta {
class Vector : public std::vector<double> {
 public:
  Vector() = default;
  Vector(const Vector &) = default;
  Vector &operator=(const Vector &) = default;
  explicit Vector(int size) : std::vector<double>(size){};
  explicit Vector(std::initializer_list<double> value)
      : std::vector<double>(value){};

  Vector &operator-=(Vector &&);
  Vector &operator/=(double);
  Vector operator*(double);
  friend void operator>>(std::fstream &, Vector &);
};
}  // namespace sfleta