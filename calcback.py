import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CSVFILE = "head300nmSiO2.csv"  # head = only 10 rows of data
phi_i = 70 * np.pi / 180  # converting incident angle from deg (first number) to rad
d_L = 300  # thickness of layer in nm
n_air = 1  # refractive index of air
rerange = 5  # upper limit for real part
imrange = 1  # upper limit for imaginary part
i = 0  # only look at one wavelength (row in csv)

csv = np.loadtxt(CSVFILE, usecols=(0,1,2,5,6),  delimiter=",", skiprows=1)

lsp_re = np.linspace(1, rerange, 1001)
lsp_im = np.linspace(0.01, imrange, 1001)
re, im = np.meshgrid (lsp_re, lsp_im, copy=False)
n_L = 1j * np.round(im,6) + np.round(re,6)
n_L = n_L.flatten() # create onedimensional array

def snell(phi, n1, n2):
    phi_ref = np.arcsin((np.sin(phi) * n1) / n2)
    return phi_ref

def fresnel(n1, phi1, n2, phi2):
    """Takes refractive indices and angles of two layers to calculate the amplitude reflection coefficients"""
    rs = (n1 * np.cos(phi1) - n2 * np.cos(phi2)) / (n1 * np.cos(phi1) + n2 * np.cos(phi2))
    rp = (n2 * np.cos(phi1) - n1 * np.cos(phi2)) / (n2 * np.cos(phi1) + n1 * np.cos(phi2))
    return rs, rp

def calc_rho(rs_al, rp_al, rs_ls, rp_ls, d_L, n_L, lambda_vac):
    beta = 2 * np.pi * d_L * n_L * np.cos(phi_L) / lambda_vac
    rp_L = (rp_al + rp_ls * np.exp(-2*1j*beta)) / (1 + rp_al * rp_ls * np.exp(-2 * 1j * beta))
    rs_L = (rs_al + rs_ls * np.exp(-2*1j*beta)) / (1 + rs_al * rs_ls * np.exp(-2 * 1j * beta))
    rho = rp_L / rs_L
    return rho

lambda_vac = csv[i, 0]
n_S = (csv[i, 3] + 1j * csv[i, 4])

phi_L = snell(phi_i, n_air, n_L)
phi_S = snell(phi_L, n_L, n_S)
# Fresnel equations:
# air/layer:
rs_al, rp_al = fresnel(n_air, phi_i, n_L, phi_L)
# layer/substrate:
rs_ls, rp_ls = fresnel(n_L, phi_L, n_S, phi_S)

rho_L = calc_rho(rs_al, rp_al, rs_ls, rp_ls, d_L, n_L, lambda_vac)

# psi is in our csv-file at index 1, delta at index 2 at row "i" for lambda
psi = csv[i][1]
delta = csv[i][2]
rho = np.tan(psi) * np.exp(1j * delta)
diff = abs(rho - rho_L)  # magnitude of complex number
idx = np.argmin(diff)  # index of the minimum
minimum = diff[idx]
n = n_L[idx]
print("The layer has the refractive index n_L = ", n)
