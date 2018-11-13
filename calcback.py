import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CSVFILE = "100nmSiO2.csv"  # head300nmSiO2.csv = only 10 rows of data
phi_i = 70 * np.pi / 180  # converting incident angle from deg (first number) to rad
d_L = 100  # thickness of layer in nm
n_air = 1  # refractive index of air
rerange = 5  # upper limit for real part
imrange = 1  # upper limit for imaginary part
# i = 88  # only look at one wavelength (row in csv)

csv = np.loadtxt(CSVFILE, delimiter=",", skiprows=1)

lsp_re = np.linspace(1, rerange, 1001)
lsp_im = np.linspace(0.01, imrange, 1001)
# re, im = np.meshgrid (lsp_re, lsp_im, copy=False)
# n_L = 1j * np.round(im,6) + np.round(re,6)
# n_L = n_L.flatten() # create onedimensional array
n_L = lsp_re

def snell(phi, n1, n2):
  """Calculates the refractive angle, parameters are incident angle phi, refractive index of first medium n1 and of second medium n2"""
  phi_ref = np.arcsin((n1/n2)*np.sin(phi))
  return phi_ref

def fresnel(n1, phi1, n2, phi2):
    """Takes refractive indices and angles of two layers to calculate the amplitude reflection coefficients"""
    rs = (n1 * np.cos(phi1) - n2 * np.cos(phi2)) / (n1 * np.cos(phi1) + n2 * np.cos(phi2))
    rp = (n2 * np.cos(phi1) - n1 * np.cos(phi2)) / (n2 * np.cos(phi1) + n1 * np.cos(phi2))
    return rs, rp

def calc_diff(n_L, rho_giv):
    #Snell's Law:
    phi_L = snell(phi_i, n_air, n_L)
    phi_S = snell(phi_L, n_L, n_S)
    # Fresnel equations:
    # air/layer:
    rs_al, rp_al = fresnel(n_air, phi_i, n_L, phi_L)
    # layer/substrate:
    rs_ls, rp_ls = fresnel(n_L, phi_L, n_S, phi_S)

    beta = (2 * np.pi / lambda_vac) * d_L * n_L * np.cos(phi_L)
    rp_L = (rp_al + rp_ls * np.exp(-2 * 1j * beta)) / (
        1 + rp_al * rp_ls * np.exp(-2 * 1j * beta))
    rs_L = (rs_al + rs_ls * np.exp(-2 * 1j * beta)) / (
        1 + rs_al * rs_ls * np.exp(-2 * 1j * beta))
    rho_L = rp_L / rs_L
    return abs(rho_giv - rho_L), rho_L

lambda_vac = csv[i][0]
n_S = (csv[i][3] + 1j * csv[i][4])

# psi is in our csv-file at index 1, delta at index 2 at row "i" for lambda
n_array = [] 
rho_array = []
for i, row in enumerate(csv):
    lambda_vac = csv[i][0]
    psi = csv[i][1] * (np.pi/180)
    delta = csv[i][2] * (np.pi/180)
    n_S = csv[i][3] + csv[i][4] * 1j
    rho_giv = np.tan(psi) * np.exp(1j * delta)
    diff, rho_L = calc_diff(n_L, rho_giv)
    idx = np.argmin(diff)  # index of the minimum
    minimum = min(diff)
    n_array = np.append(n_array, n_L[idx])
    rho_array = np.append(rho_array, rho_L[idx])
psi_L = np.arctan(abs(rho_array)) * 180/np.pi
delta_L = np.angle(rho_array) * 180/np.pi
# print("the layer has the refractive index n_L = " , n_array)
