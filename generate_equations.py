#!/usr/bin/env python3
"""
Generate the Euler-Lagrange equations in LaTeX format for documentation.
"""
import sys
sys.path.insert(0, '/home/user/slackline')

import sympy as sp
from src.core.lagrangians import ideal

# Initialize symbols
x = sp.Symbol("x")
y = sp.Function("y")(x)
n = sp.Function("n")(x)
m = sp.Symbol("m", positive=True)
g = sp.Symbol("g", positive=True)
K = sp.Symbol("K", positive=True)

print("=== Slackline Physics: Euler-Lagrange Equations ===\n")

# Get the Lagrangian
L = ideal(x, y, n, m, g, K)

print("1. THE LAGRANGIAN")
print("-" * 60)
print("The system has two generalized coordinates:")
print("  y(x) = vertical drop at horizontal position x")
print("  n(x) = natural (unstretched) length from anchor to position x")
print()
print("Lagrangian L = L_gravity + L_elastic:")
print()

# Extract components
y_x = y.diff(x)
n_x = n.diff(x)
gravity = m*g*y*n_x
tension = K/2*(1 + y_x**2)/n_x - K*sp.sqrt(1 + y_x**2) + K/2*n_x

print("Gravitational term:")
print(f"  L_gravity = {sp.latex(gravity)}")
print()
print("Elastic tension term:")
print(f"  L_elastic = {sp.latex(tension)}")
print()
print("Full Lagrangian:")
print(f"  L = {sp.latex(L)}")
print()

# Compute Euler-Lagrange equations
print("\n2. EULER-LAGRANGE EQUATIONS")
print("-" * 60)
print("The equations of motion are given by:")
print("  ∂L/∂y - d/dx(∂L/∂y') = 0")
print("  ∂L/∂n - d/dx(∂L/∂n') = 0")
print()

els = sp.calculus.euler.euler_equations(L, [y, n], [x])

print("For y(x):")
print(f"  {sp.latex(els[0])} = 0")
print()
print("For n(x):")
print(f"  {sp.latex(els[1])} = 0")
print()

# Solve for second derivatives
print("\n3. SECOND-ORDER FORM")
print("-" * 60)
print("Solving for second derivatives y'' and n'':")
print()

sol = sp.solve(els, [y.diff(x, 2), n.diff(x, 2)])
y_xx = sol[y.diff(x, 2)]
n_xx = sol[n.diff(x, 2)]

print("y''(x) =")
print(f"  {sp.latex(y_xx)}")
print()
print("n''(x) =")
print(f"  {sp.latex(n_xx)}")
print()

# Convert to first-order system
print("\n4. FIRST-ORDER ODE SYSTEM")
print("-" * 60)
print("Introducing a = dy/dx and b = dn/dx, we get a first-order system:")
print()
print("  dy/dx = a")
print("  dn/dx = b")
print(f"  da/dx = {sp.latex(y_xx.subs({y.diff(x): sp.Symbol('a'), n.diff(x): sp.Symbol('b')}))}")
print(f"  db/dx = {sp.latex(n_xx.subs({y.diff(x): sp.Symbol('a'), n.diff(x): sp.Symbol('b')}))}")
print()

# Mass boundary conditions
print("\n5. POINT MASS BOUNDARY CONDITIONS")
print("-" * 60)
print("At a point mass M, the derivatives have discontinuities.")
print("Using conservation of canonical momentum:")
print()

dL_dy_x = sp.diff(L, y.diff(x))
dL_dn_x = sp.diff(L, n.diff(x))

print("Canonical momentum conjugate to y:")
print(f"  p_y = ∂L/∂y' = {sp.latex(dL_dy_x)}")
print()
print("Canonical momentum conjugate to n:")
print(f"  p_n = ∂L/∂n' = {sp.latex(dL_dn_x)}")
print()
print("Jump conditions across mass M:")
print("  [p_y]_right - [p_y]_left = -Mg  (weight of slackliner)")
print("  [p_n]_right - [p_n]_left = 0     (natural length continuous)")
print()

print("\n" + "=" * 60)
print("These equations are solved numerically using scipy.integrate.RK45")
print("=" * 60)
