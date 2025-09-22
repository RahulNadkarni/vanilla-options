import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("/Users/rahulsmacbook/Documents/GitHub/vanilla-options/build/greeks.csv")

# Print column names to verify
print("Available columns:")
print(df.columns.tolist())
print("\nFirst few rows:")
print(df.head())

# Create the plots
plt.figure(figsize=(15, 10))

# Plot 1: Delta
plt.subplot(2, 3, 1)
plt.plot(df["S"], df["Delta_MC"], 'b-', label='Monte Carlo', linewidth=2, alpha=0.7)
plt.plot(df["S"], df["Delta_BS"], 'r--', label='Black-Scholes', linewidth=2)
plt.title("Delta vs Spot Price")
plt.xlabel("Spot Price (S)")
plt.ylabel("Delta")
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 2: Gamma
plt.subplot(2, 3, 2)
plt.plot(df["S"], df["Gamma_MC"], 'b-', label='Monte Carlo', linewidth=2, alpha=0.7)
plt.plot(df["S"], df["Gamma_BS"], 'r--', label='Black-Scholes', linewidth=2)
plt.title("Gamma vs Spot Price")
plt.xlabel("Spot Price (S)")
plt.ylabel("Gamma")
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 3: Vega
plt.subplot(2, 3, 3)
plt.plot(df["S"], df["Vega_MC"], 'b-', label='Monte Carlo', linewidth=2, alpha=0.7)
plt.plot(df["S"], df["Vega_BS"], 'r--', label='Black-Scholes', linewidth=2)
plt.title("Vega vs Spot Price")
plt.xlabel("Spot Price (S)")
plt.ylabel("Vega")
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 4: Rho
plt.subplot(2, 3, 4)
plt.plot(df["S"], df["Rho_MC"], 'b-', label='Monte Carlo', linewidth=2, alpha=0.7)
plt.plot(df["S"], df["Rho_BS"], 'r--', label='Black-Scholes', linewidth=2)
plt.title("Rho vs Spot Price")
plt.xlabel("Spot Price (S)")
plt.ylabel("Rho")
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 5: Theta
plt.subplot(2, 3, 5)
plt.plot(df["S"], df["Theta_MC"], 'b-', label='Monte Carlo', linewidth=2, alpha=0.7)
plt.plot(df["S"], df["Theta_BS"], 'r--', label='Black-Scholes', linewidth=2)
plt.title("Theta vs Spot Price")
plt.xlabel("Spot Price (S)")
plt.ylabel("Theta")
plt.legend()
plt.grid(True, alpha=0.3)

# Add a sixth subplot for error analysis
plt.subplot(2, 3, 6)
delta_error = abs(df["Delta_MC"] - df["Delta_BS"])
gamma_error = abs(df["Gamma_MC"] - df["Gamma_BS"])
plt.plot(df["S"], delta_error, 'g-', label='|Delta Error|', linewidth=2)
plt.plot(df["S"], gamma_error, 'm-', label='|Gamma Error|', linewidth=2)
plt.title("Absolute Errors (MC vs BS)")
plt.xlabel("Spot Price (S)")
plt.ylabel("Absolute Error")
plt.legend()
plt.grid(True, alpha=0.3)
plt.yscale('log') 

plt.tight_layout()
plt.show()

print("\nStatistics:")
print(f"Max Delta Error: {max(abs(df['Delta_MC'] - df['Delta_BS'])):.6f}")
print(f"Max Gamma Error: {max(abs(df['Gamma_MC'] - df['Gamma_BS'])):.6f}")
print(f"Max Vega Error: {max(abs(df['Vega_MC'] - df['Vega_BS'])):.6f}")
print(f"Max Rho Error: {max(abs(df['Rho_MC'] - df['Rho_BS'])):.6f}")
print(f"Max Theta Error: {max(abs(df['Theta_MC'] - df['Theta_BS'])):.6f}")