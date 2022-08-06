from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Use latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Import netCDF file
ncfile = './data/sho_newmark_beta.nc'
data = Dataset(ncfile)
var = data.variables

# Prepare Data to Plot
t = var['t'][:]
x_0 = var['x_0'][:]
v_0 = var['v_0'][:]  
a_0 = var['a_0'][:]
x_1 = var['x_1'][:]
v_1 = var['v_1'][:]
a_1 = var['a_1'][:]
x_2 = var['x_2'][:]
v_2 = var['v_2'][:]
a_2 = var['a_2'][:]

x_vec = [x_0, x_1, x_2]
v_vec = [v_0, v_1, v_2]
a_vec = [a_0, a_1, a_2]

zeta_vec = [0, 0.01, 0.02]

for i in range(3):
    zeta = zeta_vec[i]
    x = x_vec[i]
    v = v_vec[i]
    a = a_vec[i]
    E = 0.5 * v**2 + 0.5 * 200 * x**2

    # Prepare Plot
    plt.figure(figsize=(10,6), dpi=300)
    plt.title(r"Damped Harmonic Oscillator ($\zeta = " + str(zeta) + r"$)", fontsize=16)
    plt.xlabel(r'$t$', fontsize=14)
    plt.ylabel(r'$x$', fontsize=14)
    
    # Plot with Legends
    plt.plot(t, x, label=r'Position')
    
    # Other options
    plt.legend(fontsize=12)
    plt.grid()
    plt.savefig(f"plot/{i}0_position_zeta={zeta}.png", dpi=300)

    # Prepare Plot
    plt.figure(figsize=(10,6), dpi=300)
    plt.title(r"Damped Harmonic Oscillator ($\zeta = " + str(zeta) + r"$)", fontsize=16)
    plt.xlabel(r'$t$', fontsize=14)
    plt.ylabel(r'$v$', fontsize=14)
    
    # Plot with Legends
    plt.plot(t, v, label=r'Velocity')
    
    # Other options
    plt.legend(fontsize=12)
    plt.grid()
    plt.savefig(f"plot/{i}1_velocity_zeta={zeta}.png", dpi=300)

    # Prepare Plot
    plt.figure(figsize=(10,6), dpi=300)
    plt.title(r"Damped Harmonic Oscillator ($\zeta = " + str(zeta) + r"$)", fontsize=16)
    plt.xlabel(r'$t$', fontsize=14)
    plt.ylabel(r'$a$', fontsize=14)
    
    # Plot with Legends
    plt.plot(t, a, label=r'Acceleration')
    
    # Other options
    plt.legend(fontsize=12)
    plt.grid()
    plt.savefig(f"plot/{i}2_acceleration_zeta={zeta}.png", dpi=300)

    # Prepare Plot
    plt.figure(figsize=(10,6), dpi=300)
    plt.title(r"Damped Harmonic Oscillator ($\zeta = " + str(zeta) + r"$)", fontsize=16)
    plt.xlabel(r'$t$', fontsize=14)
    plt.ylabel(r'$E$', fontsize=14)
    
    # Plot with Legends
    plt.plot(t, E, label=r'Energy')
    
    # Other options
    plt.legend(fontsize=12)
    plt.grid()
    plt.savefig(f"plot/{i}3_energy_zeta={zeta}.png", dpi=300)
