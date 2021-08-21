# Numerical methods for PDEs

This repository contains implementations in C++ to solve partial differential equations.

<img src="advc.gif" alt="advc" align="center" style="zoom:60%;" />

### 1. Finite differences

Calculates first and second derivatives using centered finite differences.

* `derive.cpp` is the main code
* `derive.inp` is the input file with main parameters (number of steps $N$ and step size $DX$).

To compile and run the code:

```
$ g++ derive-inc.cpp -o derive.x
$ ./derive.x
```
Output:

`func.dat`, `der.dat`, `der2.dat` are the files with the function, first derivative and second derivative, respectively.



### 2. First order differential equation (eg, decay equation)

The code performs a numerical integration of a first order differential equation $\frac{dq}{dt}=-kq$ using three different schemes Euler, Heun and Leapfrog.



### 3. Diffusion equation

The code integrates the difussion equation $\frac{\partial \varphi}{\partial t}(x,t)=\alpha \frac{\partial^2 \varphi}{\partial x^2}(x,t)$ using the FTCS scheme (forward-time central-space), ie, temporal Euler and centered space scheme. 



### 4. Advection equation 

The code performs a numerical integration of the linear advection equation in one dimention $\frac{\partial u}{\partial t} + v \frac{\partial u}{\partial x} = 0$ using a centered space scheme and four temporal schemes Euler (FTCS), Heun, Lax and Leapfrog.

