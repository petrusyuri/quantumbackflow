# Integrable Defect and Quantum Backflow

## Introduction

This is a ongoing project, expanding [Bostelmann, Cadamuro and Lechner (2017)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.012112). 

## Functions supported

Until now, the following functions are supported:

1. FUN_TYPE = 1 -> Gaussian Function
2. FUN_TYPE = 2 -> Lorentzian Function
3. FUN_TYPE = 3 -> Rectangular Function

## Instructions

You can clone this project and compile it using gfortran.

```
$ git clone https://github.com/petrusyuri/quantumbackflow.git
$ cd quantumbackflow
$ make
$ ./backflow
```

## To-do simulations

### Parameters

```
N        = 2000
P_CUTOFF = 500
ALPHA    = (...)
X_LOW    = -1
X_HIGH   = +1
X_STEPS  = 100
FUN_TYPE = (...)
```

### Alpha Values

For each FUN_TYPE, ALPHA will be set to:

```
     0.00
+/-  0.05
+/-  0.40
+/-  1.00
+/- 50.00
```
