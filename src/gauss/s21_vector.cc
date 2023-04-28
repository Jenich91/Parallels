#include "sfleta_vector.h"

namespace sfleta {
Vector Vector::operator*(double value) {
  Vector res(*this);
  std::for_each(res.begin(), res.end(), [&](double& elem) { elem *= value; });
  return res;
}

Vector& Vector::operator/=(double value) {
  std::for_each(begin(), end(), [&](double& elem) { elem /= value; });
  return *this;
}

Vector& Vector::operator-=(Vector&& other) {
  for (auto it_l = begin(), it_r = other.begin(); it_l < end();
       ++it_l, ++it_r) {
    *it_l -= *it_r;
  }
  return *this;
}

void operator>>(std::fstream& fs, Vector& vec) {
  std::for_each(vec.begin(), vec.end(), [&](double& elem) { fs >> elem; });
}
}  // namespace sfleta