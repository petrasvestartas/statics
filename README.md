# Statics

C++ header-only library inspired by the book **Engineering Mechanics STATICS** by R. C. Hibbeler.

_This repository is created for self-learning purposes, driven by the frustration of not knowing even the basics of Statics...._

## Compilation
**Windows**

```bash
cd C:\brg\2_code\statics
mkdir build
cd build
cmake -G "Visual Studio 17 2022" ..
cmake --build . --config Release
.\Release\MyProject.exe
```

**Mac, Ubuntu** 
```bash
cd /path/to/brg/2_code/statics
mkdir build
cd build
cmake ..
cmake --build . --config Release
./MyProject
```

or 

```bash
chmod +x /home/petras/brg/2_code/statics/build.sh
sudo apt install cmake
/home/petras/brg/2_code/statics/build.sh
```

## Inspiration
<img width="783" alt="image" src="https://github.com/petrasvestartas/statics/assets/18013985/bde2fff1-ce70-49c9-9abf-495edb5b7266">

## Chapter 2 - Force Vectors

<p align="center">
  <img alt="Page 1" src="images/chapter2/Page 1.png">
</p>




<p align="center">
  <img alt="Page 1" src="images/chapter2/Page 2.png">
</p>


The example code can be found here: [here](/examples/chapter2/1_vector.cpp)

<p align="center">
  <img alt="Page 1" src="images/chapter2/Page 3.png">
</p>

```cpp
#include "core.hpp"

void main() {
    geo::Vector v(1, 1, 1);
    v.scale(2);
    v.rescale(-1.5); // unitize and scale the vector
    v.rescale(-0.5);
}
```

<p align="center">
  <img alt="Page 1" src="images/chapter2/Page 4.png">
</p>

```cpp
#include "core.hpp"

void main() {
    geo::Vector a(1, 0, 0);
    geo::Vector b(1.5, 0, 0);
    geo::Vector r = a + b;

    geo::Vector f1(-2, 2.5, 0);
    geo::Vector f2(4, 0, 0);
    geo::Vector r = f1 + f2;
}
```

<p align="center">
  <img alt="Page 1" src="images/chapter2/Page 5.png">
</p>


```cpp
#include "core.hpp"

void main() {

    // cosine law
    double f1 = 100;
    double f2 = 150;
    double theta_r = 115;
    double r = geo::Vector::cosine_law(f1, f2, theta_r);

    // sine law - angle
    double triangle_edge_length_a = 212.55;
    double angle_in_degrees_in_front_of_a = 115; 
    double triangle_edge_length_b = 150;
    double angle_in_degrees_in_front_of_b = geo::Vector::sine_law_angle(triangle_edge_length_a, angle_in_degrees_in_front_of_a, triangle_edge_length_b);

    // sine law - length
    double triangle_edge_length_a = 212.55;
    double angle_in_degrees_in_front_of_a = 115;
    double angle_in_degrees_in_front_of_b = 39.761714;
    double triangle_edge_length_b = geo::Vector::sine_law_length(triangle_edge_length_a, angle_in_degrees_in_front_of_a, angle_in_degrees_in_front_of_b);

}
```
