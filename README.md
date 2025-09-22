# Vanilla Options Pricing in a Black–Scholes World

## Overview
This project implements vanilla option pricing under the Black–Scholes model in C++.

It contains both:
- Closed-form formulas (analytical Black–Scholes prices)  
- Monte Carlo simulation engines (direct terminal simulation and Euler discretization)
- **Option Greeks** (Delta, Gamma, Vega, Rho, Theta) via analytical formulas and Monte Carlo methods

It also includes:
- Consistency tests to verify mathematical properties of the Black-Scholes model  
- Monte Carlo comparisons to show the convergence between simulation and formulas  

The project is modular in that option types, pricing methods, random number generators, and Greek calculators can be extended without modifying the core engine.


## Mathematics

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

### 3. Option Greeks

Greeks measure the sensitivity of option prices to various parameters. They are crucial for risk management and hedging.

#### Delta (Δ)
**Significance**: Measures price sensitivity to underlying asset price changes.  
**Hedging**: Number of shares needed for delta-neutral hedging.

![eq](https://latex.codecogs.com/png.latex?%5CDelta_%7Bcall%7D%20%3D%20%5Cfrac%7B%5Cpartial%20C%7D%7B%5Cpartial%20S%7D%20%3D%20e%5E%7B-dT%7DN%28d_1%29)

#### Gamma (Γ)
**Significance**: Measures Delta's sensitivity to price changes (convexity).  
**Hedging**: Gamma hedging protects against large price movements.

![eq](https://latex.codecogs.com/png.latex?%5CGamma_%7Bcall%7D%20%3D%20%5Cfrac%7B%5Cpartial%5E2%20C%7D%7B%5Cpartial%20S%5E2%7D%20%3D%20%5Cfrac%7Be%5E%7B-dT%7D%5Cphi%28d_1%29%7D%7BS%5Csigma%5Csqrt%7BT%7D%7D)

#### Vega (ν)
**Significance**: Sensitivity to volatility changes.  
**Hedging**: Measures volatility risk exposure.

![eq](https://latex.codecogs.com/png.latex?%5Cnu_%7Bcall%7D%20%3D%20%5Cfrac%7B%5Cpartial%20C%7D%7B%5Cpartial%20%5Csigma%7D%20%3D%20S%20e%5E%7B-dT%7D%20%5Cphi%28d_1%29%20%5Csqrt%7BT%7D)

#### Theta (Θ)
**Significance**: Time decay - price sensitivity as time passes.  
**Trading**: Important for options selling strategies.

![eq](https://latex.codecogs.com/png.latex?%5CTheta_%7Bcall%7D%20%3D%20-%5Cfrac%7B%5Cpartial%20C%7D%7B%5Cpartial%20T%7D%20%3D%20-%5Cfrac%7BS%20e%5E%7B-dT%7D%20%5Cphi%28d_1%29%20%5Csigma%7D%7B2%5Csqrt%7BT%7D%7D%20-%20rK%20e%5E%7B-rT%7D%20N%28d_2%29%20%2B%20dS%20e%5E%7B-dT%7D%20N%28d_1%29)

#### Rho (ρ)
**Significance**: Sensitivity to interest rate changes.  
**Importance**: More relevant for long-dated options.

![eq](https://latex.codecogs.com/png.latex?%5Crho_%7Bcall%7D%20%3D%20%5Cfrac%7B%5Cpartial%20C%7D%7B%5Cpartial%20r%7D%20%3D%20K%20T%20e%5E%7B-rT%7D%20N%28d_2%29)

---

### 4. Monte Carlo Greeks

Three methods are implemented for computing Greeks via Monte Carlo:

#### Finite Difference Method
- **Approach**: Bump parameter and reprice
- **Delta**: `(Price(S+ε) - Price(S)) / ε`
- **Gamma**: `(Price(S+ε) - 2×Price(S) + Price(S-ε)) / ε²`
- **Advantage**: Simple, works with any payoff
- **Disadvantage**: Computational cost, choice of ε

#### Pathwise Derivative Method
- **Approach**: Differentiate payoff along each path
- **Delta**: `∂Payoff/∂S = (ST/S) × I{ST > K}` for calls
- **Advantage**: Unbiased, lower variance
- **Disadvantage**: Requires continuous payoffs

#### Likelihood Ratio Method
- **Approach**: Use probability density derivatives (score function)
- **Delta**: `Payoff × (Z/(Sσ√T))`
- **Advantage**: Works with discontinuous payoffs
- **Disadvantage**: Higher variance

---

### 5. Monte Carlo Pricing

Monte Carlo approximates expectations:

![eq](https://latex.codecogs.com/png.latex?%5Ctext%7BPrice%7D%20%5Capprox%20e%5E%7B-rT%7D%20%5Cfrac%7B1%7D%7BN%7D%20%5Csum_%7Bi%3D1%7D%5EN%20%5Ctext%7BPayoff%7D%28S_T%5E%7B%28i%29%7D%29)

Steps:  
1. Generate *N* simulated stock prices ![eq](https://latex.codecogs.com/png.latex?S_T%5E%7B%28i%29%7D).  
2. Compute payoff for each path.  
3. Discount and average.  
4. Variance gives error bars (standard error ~ ![eq](https://latex.codecogs.com/png.latex?1/%5Csqrt%7BN%7D)).  

---

### 6. Euler Discretization
Instead of simulating directly to maturity, split into steps of size Δt:

![eq](https://latex.codecogs.com/png.latex?S_%7Bt+%5CDelta%20t%7D%20%3D%20S_t%20%2B%20%28r-d%29S_t%5CDelta%20t%20%2B%20%5Csigma%20S_t%20%5Csqrt%7B%5CDelta%20t%7D%20Z)

This is essential for **path-dependent options** (e.g., Asian options).

---

## Implementation


### Key Components
- `BlackScholes.*` → closed-form formulas  
- `Options.*` → option payoff classes  
- `MonteCarlo.*` → terminal simulation engine + Greeks calculators
- `EulerMonteCarlo.*` → pathwise simulation engine  
- `Greeks.*` → analytical Greeks formulas  # NEW
- `Utils.*` → random number generation  
- `main.cpp` → demo program  
- `tests/` → mathematical consistency + simulation vs formula convergence  

### Greek Calculation Methods
The project implements multiple approaches:

**Analytical Methods** (in `Greeks.cpp`):
- Exact Black-Scholes formulas for maximum accuracy
- Used as benchmarks for Monte Carlo methods

**Monte Carlo Methods** (in `MonteCarlo.cpp`):
- Finite Differences: General-purpose, bump-and-reprice
- Pathwise Derivatives: Efficient for smooth payoffs  
- Likelihood Ratio: Robust for discontinuous payoffs

---

## Running the Project

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


## Applications 

### Risk Management 

- Delta Hedging: Create delta-neutral portfolios
- Gamma Scalping: Profit from large price movements
- Vega Trading: Take positions on volatility changes

## Trading Strategies 

- Calendar Spreads: Theta decay exploitation
- Volatility Arbitrage: Vega mispricing opportunities
- Dynamic Hedging: Continuous portfolio rebalancing