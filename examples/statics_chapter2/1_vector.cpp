#include "core.hpp"

int main() {
  geo::Vector v(0.57, -158.63, 180.890);
  double x = v[0];
  double y = v[1];
  double z = v[2];
  double length = v.length();
  geo::log(v.to_string());
  return 1;
}