import numpy as np
import matplotlib.pyplot as plt

# Constants
E = 511e3  # Photon energy in eV
m_ec2 = 511e3
re = 2.818e-15

# Î¸ range
theta_deg = np.linspace(0, 180, 1000)
theta_rad = np.deg2rad(theta_deg)

# Scattered photon energy E'
E_prime = E / (1 + (E / m_ec2) * (1 - np.cos(theta_rad)))

# Ï = 0Â°: assume polarization aligned with scattering plane (if unpolarized, use sinÂ²Î¸ term only)
phi_rad = 0
ratio = E_prime / E
prefactor = (re**2) / 2

# Polarized Klein-Nishina differential cross-section
dsigma_domega = prefactor * ratio**2 * (
    ratio + 1/ratio - 2 * np.sin(theta_rad)**2 * np.cos(phi_rad)**2
)

# Normalize to get a probability density
pdf = dsigma_domega / np.trapz(dsigma_domega * np.sin(theta_rad), theta_rad)

# Plot
plt.figure(figsize=(8, 6))
plt.plot(theta_deg, pdf, label='Normalized PDF (single scatter)', color='purple')
plt.xlabel('Scattering Angle Î¸ (degrees)')
plt.ylabel('Probability Density (1/sr)')
plt.title('Angular Distribution of Single Compton Scattering (511 keV)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
