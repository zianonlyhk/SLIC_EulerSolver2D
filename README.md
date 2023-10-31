# Solving the 2D Euler equations using the SLIC scheme

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Usage](#usage)

## About <a name = "about"></a>

This is one of the CCM1/ACM1 practical assignments during my MPhil Scientific Computing course at the University of Cambridge. It marks my progress learning C++ and other tools in the developer suite.

## Getting Started <a name = "getting_started"></a>

Simply clone the repo and put it somewhere in your computer. Things are portable.

### Prerequisites

C++ compiler and Make software are required but they are shipped by default in most of the OS out there. The only thing special is `gnuplot` if one wants to use the codes to the fullest extent.

## Usage <a name = "usage"></a>

I first edit the absolute directory of the repo on line 128 of the file "`src/exec_interface.cc`" with the correct value, in this case it's "`/Users/zian/Room214N/github/SLIC_EulerSolver2D/`". Notice that it ends with a "`/`".

I will then type the following command on my mac terminal under the repo directory:

````
make && ./bin/slic_solver_2d && ./gif_plot.gp
````

An example set of density solutions is presented as:
![](https://raw.githubusercontent.com/zianonlyhk/SLIC_EulerSolver2D/main/rhoPlot.gif)
