# Numerical methods for PDEs

This repository contains implementations in C++ to solve partial differential equations.

| ![advc](advc.gif) |
|:--:|
| *Linear Advection equation in 2D using Euler time integration* |


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

The code performs a numerical integration of a first order differential equation  
![\Large \frac{dq}{dt}=-k](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bdq%7D%7Bdt%7D%3D-k)  
using three different schemes Euler, Heun and Leapfrog.
To compile and run the code:

```
$ g++ temp-scheme.cpp -o temp-scheme.x
$ ./temp-scheme.x
```

**Second order**

The code performs a numerical integration of a second order differential equation  
![\Large \frac{d^2q}{dt^2}=-\omega_0^2q](https://latex.codecogs.com/svg.latex?%5Cfrac%7Bd%5E2q%7D%7Bdt%5E2%7D%3D-%5Comega_0%5E2q)  

For each integration step it calculates the system energy using  
![\Large E=\frac{1}{2}m\dot{q}^2+\frac{1}{2}m\omega_0^2q^2](https://latex.codecogs.com/svg.latex?E%3D%5Cfrac%7B1%7D%7B2%7Dm%5Cdot%7Bq%7D%5E2&plus;%5Cfrac%7B1%7D%7B2%7Dm%5Comega_0%5E2q%5E2).  
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

The code integrates the difussion equation  
![\Large \frac{\partial \varphi}{\partial t}(x,t)=\alpha \frac{\partial^2 \varphi}{\partial x^2}(x,t)](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cpartial%20%5Cvarphi%7D%7B%5Cpartial%20t%7D%28x%2Ct%29%3D%5Calpha%20%5Cfrac%7B%5Cpartial%5E2%20%5Cvarphi%7D%7B%5Cpartial%20x%5E2%7D%28x%2Ct%29)  
using the FTCS scheme (forward-time central-space), ie, temporal Euler and centered space scheme. 

To compile and run the code:

```
$ g++ FTCS.cpp -o FTCS.x
$ ./FTCS.x
```

### 4. Advection equation

**1D advection equation**

The code performs a numerical integration of the linear advection equation in one dimention  
![\Large \frac{\partial u}{\partial t} + v \frac{\partial u}{\partial x} = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20t%7D%20&plus;%20v%20%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20x%7D%20%3D%200)  
using a centered space scheme and four temporal schemes Euler (FTCS), Heun, Lax and Leapfrog.

To compile and run the code:

```
$ g++ adv-1D.cpp -o adv-1D.x
$ ./adv-1D.x
```

**2D advection equation**

The code performs a numerical integration of the 2D advection equation  
![\Large \frac{\partial u}{\partial t}(x,y,t) + \vec{v} \vec{\nabla} u(x,y,t) = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20t%7D%28x%2Cy%2Ct%29%20&plus;%20%5Cvec%7Bv%7D%20%5Cvec%7B%5Cnabla%7D%20u%28x%2Cy%2Ct%29%20%3D%200)   


![\Large \frac{\partial u}{\partial t} + v_x \frac{\partial u}{\partial x} + v_y \frac{\partial u}{\partial y} = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20t%7D%20&plus;%20v_x%20%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20x%7D%20&plus;%20v_y%20%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20y%7D%20%3D%200)  
using an FTCS scheme.

To compile and run the code:

```
$ g++ adv-2D.cpp -o adv-2D.x
$ ./adv-1D.x
```