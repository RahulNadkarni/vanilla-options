# Vanilla Options Pricing in a Black–Scholes World

## Overview
This project implements vanilla option pricing under the Black–Scholes model in C++.

It contains both:
- Closed-form formulas (analytical Black–Scholes prices)  
- Monte Carlo simulation engines (direct terminal simulation and Euler discretization)

It also includes:
- Consistency tests to verify mathematical properties of the Black-Scholes model  
- Monte Carlo comparisons to demonstrate convergence between simulation and formulas  

The project is modular — option types, pricing methods, and random number generators can be extended without modifying the core engine.

---

## The Math Behind It

### 1. Stock Dynamics
In a risk-neutral world, the stock follows a geometric Brownian motion:

![eq](https://latex.codecogs.com/png.latex?dS_t%20%3D%20%28r-d%29S_tdt%20%2B%20%5Csigma%20S_tdW_t)

Where:  
- *r* = risk-free rate  
- *d* = continuous dividend yield  
- *σ* = volatility  
- *W<sub>t</sub>* = Brownian motion  

Closed-form solution for terminal price:

![eq](https://latex.codecogs.com/png.latex?S_T%20%3D%20S_0%20%5Cexp%5Cleft%28%28r-d-%5Cfrac%7B1%7D%7B2%7D%5Csigma%5E2%29T%20%2B%20%5Csigma%20%5Csqrt%7BT%7D%20Z%5Cright%29%2C%20%5Cquad%20Z%20%5Csim%20N%280%2C1%29)

---

### 2. Black–Scholes Formulas

**PDF of standard normal**  
![eq](https://latex.codecogs.com/png.latex?%5Cphi%28x%29%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%7D%7D%20e%5E%7B-%5Cfrac%7B1%7D%7B2%7Dx%5E2%7D)

**CDF of standard normal**  
![eq](https://latex.codecogs.com/png.latex?N%28x%29%20%3D%20%5Cint_%7B-%5Cinfty%7D%5Ex%20%5Cphi%28t%29dt)

**Intermediate terms**  

![eq](https://latex.codecogs.com/png.latex?d_%7B1%2C2%7D%20%3D%20%5Cfrac%7B%5Cln%28S_0%2FK%29%20%2B%20%28r-d%20%5Cpm%20%5Cfrac%7B1%7D%7B2%7D%5Csigma%5E2%29T%7D%7B%5Csigma%5Csqrt%7BT%7D%7D)

**Formulas:**  

- Forward price  
![eq](https://latex.codecogs.com/png.latex?F%20%3D%20S_0e%5E%7B-dT%7D%20-%20Ke%5E%7B-rT%7D)

- Call  
![eq](https://latex.codecogs.com/png.latex?C%20%3D%20S_0e%5E%7B-dT%7DN%28d_1%29%20-%20Ke%5E%7B-rT%7DN%28d_2%29)

- Put  
![eq](https://latex.codecogs.com/png.latex?P%20%3D%20Ke%5E%7B-rT%7DN%28-d_2%29%20-%20S_0e%5E%7B-dT%7DN%28-d_1%29)

- Digital Call  
![eq](https://latex.codecogs.com/png.latex?DC%20%3D%20e%5E%7B-rT%7DN%28d_2%29)

- Digital Put  
![eq](https://latex.codecogs.com/png.latex?DP%20%3D%20e%5E%7B-rT%7DN%28-d_2%29)

- Zero-Coupon Bond  
![eq](https://latex.codecogs.com/png.latex?B%20%3D%20e%5E%7B-rT%7D)

---

### 3. Monte Carlo Pricing

Monte Carlo approximates expectations:

![eq](https://latex.codecogs.com/png.latex?%5Ctext%7BPrice%7D%20%5Capprox%20e%5E%7B-rT%7D%20%5Cfrac%7B1%7D%7BN%7D%20%5Csum_%7Bi%3D1%7D%5EN%20%5Ctext%7BPayoff%7D%28S_T%5E%7B%28i%29%7D%29)

Steps:  
1. Generate *N* simulated stock prices ![eq](https://latex.codecogs.com/png.latex?S_T%5E%7B%28i%29%7D).  
2. Compute payoff for each path.  
3. Discount and average.  
4. Variance gives error bars (standard error ~ ![eq](https://latex.codecogs.com/png.latex?1/%5Csqrt%7BN%7D)).  

---

### 4. Euler Discretization
Instead of simulating directly to maturity, split into steps of size Δt:

![eq](https://latex.codecogs.com/png.latex?S_%7Bt+%5CDelta%20t%7D%20%3D%20S_t%20%2B%20%28r-d%29S_t%5CDelta%20t%20%2B%20%5Csigma%20S_t%20%5Csqrt%7B%5CDelta%20t%7D%20Z)

This is essential for **path-dependent options** (e.g., Asian options).

---

## Implementation

### Directory Structure
 project_root/
├── CMakeLists.txt
├── include/
│ ├── BlackScholes.h
│ ├── Options.h
│ ├── MonteCarlo.h
│ ├── EulerMonteCarlo.h
│ └── Utils.h
├── src/
│ ├── main.cpp
│ ├── BlackScholes.cpp
│ ├── Options.cpp
│ ├── MonteCarlo.cpp
│ ├── EulerMonteCarlo.cpp
│ └── Utils.cpp
└── tests/
├── ConsistencyTests.cpp
└── MonteCarloComparison.cpp


### Key Components
- `BlackScholes.*` → closed-form formulas  
- `Options.*` → option payoff classes  
- `MonteCarlo.*` → terminal simulation engine  
- `EulerMonteCarlo.*` → pathwise simulation engine  
- `Utils.*` → random number generation  
- `main.cpp` → demo program  
- `tests/` → mathematical consistency + simulation vs formula convergence  

---

## ▶️ Running the Project

### 1. Build
```bash
mkdir build && cd build
cmake ..
make 
```

### 2. Run Main Program 

```bash 
    ./vanilla_options
```
### 3. Run Consistency Checks 

```bash 
    ./consistency_tests
```
Checks:
Put–call parity
Call decreasing in strike
Call increasing in volatility
Convexity in strike
Digital call + put = bond




