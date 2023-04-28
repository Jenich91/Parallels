#pragma once
#include <climits>
#include <functional>
#include <map>
#include <memory>
#include <sstream>

#include "ant/sfleta_ant_algorithm.h"
#include "gauss/sfleta_gauss.h"
#include "sfleta_timer.h"
#include "winograd/sfleta_winograd_algorithm.h"

namespace sfleta {
class Interface {
 public:
  ~Interface() {}
  Interface(const Interface &) = delete;
  Interface(const Interface &&) = delete;
  Interface &operator=(const Interface &) = delete;
  Interface &operator=(const Interface &&) = delete;

  static Interface &GetInstance() {
    static Interface instance;
    return instance;
  }

  void Show(int input = 0);

 private:
  Interface();
  void WinogradFromFile(Matrix *&, Matrix *&);
  void WinogradFromRandom(Matrix *&, Matrix *&);
  void WaitingForInput();
  void WrongInputAttention();
  void PrintAntMenu();

  void AntStart();
  void GaussStart();
  void WinogradStart();

  std::map<int, std::function<void()> > dictionary;
};
}  // namespace sfleta
