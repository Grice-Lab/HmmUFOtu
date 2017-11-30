EGriceLab::util introduction
============================
EGriceLab::util is shared C++ library used among EGriceLab which contains useful classes and functions
for common functionalities across multiple C/C++ projects

Implementation
--------------
EGriceLab::util is written in pure C++98, and built with the GNU Autotools

Dependencies
------------
EGriceLab::util depends on the popular head-only C++ library Boost.

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
