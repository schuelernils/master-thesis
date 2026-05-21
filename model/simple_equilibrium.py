import cantera as ct

T = 1100 + 273.15  # K
p = 10 * 1e5       # Pa
X_feed = {'CH4': 0.25, 'H2O': 0.5, 'O2': 0.05, 'N2': 0.2}

g = ct.Solution('gri30.yaml')
g.TPX = T, p, X_feed
w_CH4_in = g['CH4'].Y[0]

g.equilibrate('TP')
w_CH4_eq = g['CH4'].Y[0]

X_CH4 = (w_CH4_in - w_CH4_eq) / w_CH4_in
print(f"CH4-Umsatz bei T={T-273.15:.0f} °C, p={p/1e5:.0f} bar: {X_CH4:.3f}")
