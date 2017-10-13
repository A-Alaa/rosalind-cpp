# rosalind-cpp
This project aims to implement bioinformatics algorithms (BA track at rosalind.info) making use of C++11/14 STL.
The Catch.hpp header-only library is used to write test cases for most algorithms. This project might be used further as an auxillary computational library in advanced projects.
Currently ba1*, ba2*, ba3*, and ba4* are completely implemented. The rest are yet to come.

## Instruction 
### Requirements 
* CMake 3.4  

### Installation
* get the repository.
```
git clone https://github.com/A-Alaa/rosalind-cpp.git
```
* Building:
```
cd rosalind-cpp
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=<Release|Debug>..
make -j4
make install
```

## Author 
Asem Alaa  
rosalind.info: http://rosalind.info/users/asem_alla/  
email: asem_alla@msn.com
