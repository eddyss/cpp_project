#include <iostream>
#include "var_alea.hpp"

int main() {
  init_alea();

  uniform U(0,1);
  std::cout << U() << std::endl;
  std::cout << U() << std::endl;
  std::cout << "JP" << std::endl;

  return 0;
};
