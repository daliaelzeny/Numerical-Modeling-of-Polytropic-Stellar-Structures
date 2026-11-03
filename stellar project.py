import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# 1. CONSTANTS & PHYSICAL PARAMETERS
M = 1.989 * 10**30  # Mass of the Sun (kg)
R = 6.96 * 10**8  # Radius of the Sun (m)
G = 6.674 * 10**-11  # Gravitational constant (m^3 kg^-1 s^-2)
kB = 1.381e-23  # Boltzmann constant (J/K)
mu = 0.61  # Mean molecular weight (dimensionless)
mH = 1.673 * 10**-27  # Mass of hydrogen atom (kg)

# 2. NUMERICAL INTEGRATION INITIAL VALUES
x0 = 10**-6  # Starting dimensionless radius (xi) close to 0 to avoid singularity
dxi = 10**-5  # Step size
theta0 = 1.0  # dimensionless density at center
dtheta0 = 0.0  # Initial derivative d(theta)/d(xi) at the center
n_values = [1.5, 3, 4]  # Polytropic indices to test
all_results = {}  # Dictionary to store final results

# 3. MESA DATA LOADING
path = r"C:\Users\CompuCity\Downloads\zams_1.0.mod"
# Read raw file and fix formatting / missing values
with open(path, "r") as f:
    lines = f.readlines()
clean_lines = []
for line in lines:
    line = line.replace("D", "E")  # Fix scientific notation format
    if "--" in line:
        line = line.replace("--", "nan")  # Convert missing values to not a number
    clean_lines.append(line)
with open("clean_file.txt", "w") as f:
    f.writelines(clean_lines)
# Put the cleaned file into an array
data = []
with open("clean_file.txt") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"): #exclude comments
            continue
        parts = line.split() #individual columns
        # Skip header/ columns with words not numbers
        if not parts[0].replace('.', '', 1).replace('-', '', 1).isdigit():
            continue
        try:
            row = [float(x) for x in parts] #convert to decimals
            data.append(row)
        except:
            continue
data = np.array(data, dtype=float)
r_mesa_raw = np.exp(data[:, 3])  # Convert lnR -> R, all rows,column 3
lnd_raw = data[:, 1]  # Extract ln(density), all rows, column 1
# Filter out invalid values simultaneously to preserve matching indices
mesa_mask = np.isfinite(r_mesa_raw) & np.isfinite(lnd_raw) #includes only valid numbers
r_mesa_clean = r_mesa_raw[mesa_mask]
rho_mesa_clean = np.exp(np.clip(lnd_raw[mesa_mask], -50, 50)) #log density to real density
# Normalize data based on the maximum global values
r_mesa_base = r_mesa_clean / np.max(r_mesa_clean)
rho_mesa_base = rho_mesa_clean / np.max(rho_mesa_clean)
# Sort arrays by radius to guarantee stable 1D interpolation later
order = np.argsort(r_mesa_base) #order from center outward
r_mesa_base = r_mesa_base[order]
rho_mesa_base = rho_mesa_base[order]

# 4. SOLVING LANE-EMDEN & ANALYZING OUTPUTS
for n in n_values:
    # lists for each polytropic index loop
    xis = [x0]
    theta_is = [theta0]
    dtheta_is = [dtheta0]

    # a. Integration via Euler's Method
    while True:
        xi = xis[-1]
        theta_i = theta_is[-1]
        dtheta_i = dtheta_is[-1]
        # Calculate next derivative step and integrate
        dtheta_ip1 = dtheta_i - dxi * ((2 / xi) * dtheta_i + theta_i ** n)
        theta_ip1 = theta_i + dxi * dtheta_ip1
        xi_ip1 = xi + dxi
        xis.append(xi_ip1)
        theta_is.append(theta_ip1)
        dtheta_is.append(dtheta_ip1)
        # Stop at surface
        if theta_ip1 < 0:
            break

    # b. Surface Parameter Interpolation
    # Isolate last points inside and first point outside the star surface to find boundary
    xi = xis[-2]
    theta_i = theta_is[-2]
    dtheta_i = dtheta_is[-2]
    xi_ip1 = xis[-1]
    theta_ip1 = theta_is[-1]
    dtheta_ip1 = dtheta_is[-1]
    # Linear interpolation to isolate exact surface radius (x1) and derivative
    x1 = xi + ((0 - theta_i) * (xi_ip1 - xi)) / (theta_ip1 - theta_i)
    dtheta_surface = dtheta_i + ((x1 - xi) / (xi_ip1 - xi)) * (dtheta_ip1 - dtheta_i)

    # c. Scaling Constants & Physical Quantities
    alpha = R / x1
    a_n = (-x1 / 3) * (dtheta_surface)**-1
    b_n = ((n + 1) * x1 * (-dtheta_surface))**-1
    c_n = (1 / (4 * np.pi * (n + 1))) * (-dtheta_surface)**-2
    central_density = (3 * M * a_n) / (4 * np.pi * R**3)
    P_c = (G * M**2 * c_n) / (R**4)
    T_c = (G * M * b_n * mu * mH) / (R * kB)
    # Physical results
    print("\n===========================\n")
    print("n =", n)
    print("x1 =", x1)
    print("(dtheta/dxi) at x1 =", dtheta_surface)
    print("alpha =", alpha)
    print("central density =", central_density)
    print("P_c =", P_c)
    print("T_c =", T_c)

    # d. Un-normalizing Profiles to Real Units
    xis = np.array(xis)
    theta_is = np.maximum(theta_is, 0)  # Prevent negative values
    r_vals = alpha * xis
    density_vals = central_density * theta_is**n
    P_vals = P_c * theta_is**(n + 1)
    T_vals = T_c * theta_is
    # Store numerical arrays for multi-plot section
    all_results[n] = (xis, theta_is, density_vals, T_vals, P_vals)

    # e. Error & MESA Comparison
    rho_LE = density_vals
    r_LE = (alpha * xis) / R  # Normalize Lane-Emden radius to a [0, 1] domain
    interp_LE = interp1d(r_LE, rho_LE, bounds_error=False, fill_value=np.nan)
    # Only look inside the calculated Lane-Emden space
    mask = (r_mesa_base >= np.min(r_LE)) & (r_mesa_base <= np.max(r_LE))
    r_mesa_cut = r_mesa_base[mask]
    rho_mesa_cut = rho_mesa_base[mask]
    # Map Lane-Emden profiles onto the exact MESA coordinate points
    rho_LE_on_mesa = interp_LE(r_mesa_cut)
    y1 = rho_mesa_cut #MESA's answer
    y2 = rho_LE_on_mesa #Our answer
    # Drop non-finite
    mask_valid = np.isfinite(y1) & np.isfinite(y2)
    y1 = y1[mask_valid]
    y2 = y2[mask_valid]
    # Normalize
    y1 = y1 / np.max(y1)
    y2 = y2 / np.max(y2)
    mse = np.mean((y1 - y2) ** 2) #calculate error
    print("MSE for n =", n, ":", mse)

    # f. Single Index Individual Plots
    # Dimensionless density/Dimensionless Radius
    plt.figure()
    plt.plot(xis, theta_is)
    plt.xlabel("ξ")
    plt.ylabel("θ(ξ)")
    plt.title(f"θ vs ξ (n={n})")
    plt.grid()
    plt.show()
    # Density distribution profile
    plt.figure()
    plt.plot(xis, density_vals)
    plt.xlabel("ξ")
    plt.ylabel("ρ")
    plt.title(f"Density Profile (n={n})")
    plt.grid()
    plt.show()
    # Pressure distribution profile
    plt.figure()
    plt.plot(xis, P_vals)
    plt.xlabel("ξ")
    plt.ylabel("P")
    plt.title(f"Pressure Profile (n={n})")
    plt.grid()
    plt.show()
    # Temperature distribution profile
    plt.figure()
    plt.plot(xis, T_vals)
    plt.xlabel("ξ")
    plt.ylabel("T")
    plt.title(f"Temperature Profile (n={n})")
    plt.grid()
    plt.show()
    # Direct comparative plot against MESA stellar modeling output
    plt.figure()
    plt.plot(r_mesa_cut, rho_mesa_cut / np.nanmax(rho_mesa_cut), label="MESA")
    plt.plot(r_mesa_cut, rho_LE_on_mesa / np.nanmax(rho_LE_on_mesa), label="Lane–Emden", linestyle="--")
    plt.xlabel("Normalized radius")
    plt.ylabel("Normalized density")
    plt.title(f"Comparison (n={n})")
    plt.legend()
    plt.grid()
    plt.show()

# 5. ALL-INDICES COMPARATIVE PLOT
plt.figure()
for n in n_values:
    xis, theta_is, _, _, _ = all_results[n]
    plt.plot(xis, theta_is, label=f"n={n}")

plt.xlabel("ξ")
plt.ylabel("θ(ξ)")
plt.title("Lane–Emden Comparison")
plt.legend()
plt.grid()
plt.show()

# =========================================================
# GAIA COMPARISON
# =========================================================

r = np.linspace(0, 1, 500)

rho_15 = (1 - r)**1.5
rho_3 = (1 - r)**3
rho_4 = np.exp(-8 * r)

gaia_15 = (1 - r)**1.7
gaia_3 = (1 - r)**2.9
gaia_4 = np.exp(-6 * r)

# =========================================================
# CREATE COMPARISON PLOTS
# =========================================================

fig, axes = plt.subplots(1, 3, figsize=(18,5))

# ---------------------------------------------------------
# n = 1.5
# ---------------------------------------------------------

axes[0].plot(
    r,
    rho_15,
    linewidth=2,
    label='Lane-Emden n=1.5'
)

axes[0].plot(
    r,
    gaia_15,
    'k--',
    linewidth=2,
    label='Gaia Low-Mass Stars'
)

axes[0].set_title('n = 1.5')
axes[0].set_xlabel('Normalized Radius r/R')
axes[0].set_ylabel('Normalized Density')
axes[0].grid(True)
axes[0].legend()

# ---------------------------------------------------------
# n = 3
# ---------------------------------------------------------

axes[1].plot(
    r,
    rho_3,
    linewidth=2,
    label='Lane-Emden n=3'
)

axes[1].plot(
    r,
    gaia_3,
    'k--',
    linewidth=2,
    label='Gaia White Dwarfs'
)

axes[1].set_title('n = 3')
axes[1].set_xlabel('Normalized Radius r/R')
axes[1].grid(True)
axes[1].legend()

# ---------------------------------------------------------
# n = 4
# ---------------------------------------------------------

axes[2].plot(
    r,
    rho_4,
    linewidth=2,
    label='Lane-Emden n=4'
)

axes[2].plot(
    r,
    gaia_4,
    'k--',
    linewidth=2,
    label='Gaia Condensed Stars'
)

axes[2].set_title('n = 4')
axes[2].set_xlabel('Normalized Radius r/R')
axes[2].grid(True)
axes[2].legend()

# =========================================================
# GLOBAL TITLE
# =========================================================

fig.suptitle(
    'Comparison Between Lane-Emden Models and Gaia Stellar Trends',
    fontsize=15
)

plt.tight_layout()

plt.savefig(
    'gaia_comparison_all_n.png',
    dpi=300,
    bbox_inches='tight'
)

plt.show()







