# Raman

Simulate 1D/2D Raman signals via the hierarchical equation of motion approach.

## Dependence

* [Lapack](http://www.netlib.org/lapack/)

* [OpenBlas](http://www.openblas.net/) (or [Blas](http://www.netlib.org/blas/))

* [Armadillo](http://arma.sourceforge.net/)

* [json11](https://github.com/dropbox/json11)

## Install

1. git clone https://github.com/hou-dao/raman.git raman

2. cd raman

3. mkdir build && cmake ..

4. make

The binary files are located in raman/bin.

## Usage

1. cd raman/input

2. update gen_input.py and run it to generate input files

3. start simulation with ../bin/nonlinear

## Contact

[Zhang Houdao](mailto:houdao@connect.ust.hk)
