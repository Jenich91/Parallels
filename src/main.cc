#include <iostream>

#include "sfleta_interface.h"

int main(int argc, const char** argv) {
  if (argc == 2) {
    std::map<std::string, int> algorithm_names;
    algorithm_names["ant"] = 1;
    algorithm_names["gauss"] = 2;
    algorithm_names["winograd"] = 3;

    std::cout << "~" << argv[1] << " started"
              << "~" << std::endl;
    sfleta::Interface::GetInstance().Show(algorithm_names[argv[1]]);
  }
  return 0;
}
