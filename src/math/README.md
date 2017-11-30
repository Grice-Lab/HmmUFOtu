EGriceLab::math introduction
============================
EGriceLab::math is shared C++ math library used among EGriceLab which contains useful classes and functions
for common functionalities across multiple C/C++ projects

Implementation
--------------
EGriceLab::math is written in pure C++98, and built with the GNU Autotools

Dependencies
------------
EGriceLab::math depends on the popular header-only C++ libraries Boost and Eigen3.
Some of its functionality is dependent on Boost::Math and thus it is optional to build/install Boost.

Installation
------------
1. Configure installation, by running the command
```bash
./configure
```
You may consider providing additional options, such as `--prefix`, `--with-boost`, etc.

2. Compile and link, by running the command
```bash
make
```

3. Install
```bash
make install
```
You may need root privilege to do it, such as using `sudo`.
