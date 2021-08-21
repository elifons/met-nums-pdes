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

### 2. ODEs

The following codes perform numerical integration of Ordinary Differential Equations. In all cases `cond.inp` is the input file with main parameters.

**First order**

The code performs a numerical integration of a first order differential equation ![\Large \frac{dq}{dt}=-k](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdq%7D%7Bdt%7D%3D-k) using three different schemes Euler, Heun and Leapfrog.
To compile and run the code:

```
$ g++ temp-scheme.cpp -o temp-scheme.x
$ ./temp-scheme.x
```

**Second order**

The code performs a numerical integration of a second order differential equation ![\Large \frac{d^2q}{dt^2}=-\omega_0^2q](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bd%5E2q%7D%7Bdt%5E2%7D%3D-%5Comega_0%5E2q). For each integration step it calculates the system energy using ![\Large E=\frac{1}{2}m\dot{q}^2+\frac{1}{2}m\omega_0^2q^2](https://latex.codecogs.com/svg.latex?E%3D%5Cfrac%7B1%7D%7B2%7Dm%5Cdot%7Bq%7D%5E2&plus;%5Cfrac%7B1%7D%7B2%7Dm%5Comega_0%5E2q%5E2).
To compile and run the code:

```
$ g++ temp-scheme-SO.cpp -o temp-scheme-SO.x
$ ./temp-scheme-SO.x
```
**Pendulum**
The code performs a numerical integration of the pendulum equation.
To compile and run the code:

```
$ g++ pendulum.cpp -o pendulum.x
$ ./pendulum.x
```


### 3. Diffusion equation

The code integrates the difussion equation $\frac{\partial \varphi}{\partial t}(x,t)=\alpha \frac{\partial^2 \varphi}{\partial x^2}(x,t)$ using the FTCS scheme (forward-time central-space), ie, temporal Euler and centered space scheme. 



### 4. Advection equation 

The code performs a numerical integration of the linear advection equation in one dimention $\frac{\partial u}{\partial t} + v \frac{\partial u}{\partial x} = 0$ using a centered space scheme and four temporal schemes Euler (FTCS), Heun, Lax and Leapfrog.

