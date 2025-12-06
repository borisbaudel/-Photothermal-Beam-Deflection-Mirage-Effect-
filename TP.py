import numpy as np
import matplotlib.pyplot as plt

# Data
y = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,
              1,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])
I = np.array([-830,-827,-827,-824,-807,-764,-682,-551,-458,-354,-219,-85,66,
              189,335,461,580,672,760,822,910,949,963,966,967,967,968])

# Select linear region manually
mask = (y >= 0.8) & (y <= 1.2)
y_lin = y[mask]
I_lin = I[mask]

# Least-squares linear fit
coef = np.polyfit(y_lin, I_lin, 1)
slope = coef[0]
intercept = coef[1]

print("Pente (sensibilité S2) =", slope)

# Plot
plt.figure(figsize=(7,5))
plt.plot(y, I, 'o', label="Données")
plt.plot(y_lin, slope*y_lin + intercept, '-', 
         label=f"Ajustement linéaire (pente = {slope:.2f})")
plt.xlabel("y (mm)")
plt.ylabel("I (Microampères)")
plt.title("Courbe I(y) avec ajustement linéaire")
plt.grid(True)
plt.legend()
plt.show()