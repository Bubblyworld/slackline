# Slackline Physics Simulator

Physics-based slackline simulator using Lagrangian mechanics to model elastic behavior under load.

## Installation

```bash
poetry install
```

## Usage

### Python API

```python
from src.api import Constraints

# Create a 50m slackline with 2000N standing tension
constraints = Constraints(gap_length=50, anchor_tension=2000)

# Add a 75kg person at the 25m mark
constraints.add_slackliner(position=25, mass=75)

# Compute the rig
rig = constraints.rig()

# Access results: x, y (curve), T (tension), A (angle), n/l (lengths)
print(f"Max drop: {min(rig.y):.2f}m, Max tension: {max(rig.T):.0f}N")
```

### Flask REST API

```bash
python src/server.py
```

```bash
curl -X POST http://localhost:5000/rig \
  -H "Content-Type: application/json" \
  -d '{"gap_length": 50, "anchor_tension": 2000, "slackliners": [[25, 75]]}'
```

## Examples

### 25m line with 80kg person at center

![25m](plot1_25m.png)

### 100m line with two people (70kg at 30m, 80kg at 70m)

![100m](plot2_100m.png)

### 500m line with 75kg person at 100m

![500m](plot3_500m.png)

## How it works

### The Physics

Slacklines are modeled as elastic continua using **Lagrangian mechanics**. The system has two generalized coordinates:
- `y(x)` - vertical drop at horizontal position `x`
- `n(x)` - natural (unstretched) length of webbing from anchor to position `x`

The arc length element is `dl = √(1 + y'²) dx`, and the stretch ratio is `dl/dn`. For Hookean elastic webbing with stiffness `K`, the strain energy density is `K/2 · (dl/dn - 1)²`.

The **Lagrangian** combines gravitational potential energy and elastic energy:

```
L = m g y n' + K/2 · (1 + y'²)/n' - K√(1 + y'²) + K/2 · n'
```

where:
- `m` = mass per meter of webbing (kg/m)
- `g` = gravitational acceleration (9.81 m/s²)
- `K` = stiffness (N per 100% strain)
- Primes denote derivatives with respect to `x`

### The Equations

Applying the Euler-Lagrange equations `∂L/∂q - d/dx(∂L/∂q') = 0` for each coordinate yields two coupled second-order ODEs. These are symbolically derived using SymPy and converted to a first-order system:

```
dy/dx = a
dn/dx = b
da/dx = f₁(x, y, n, a, b; m, g, K)
db/dx = f₂(x, y, n, a, b; m, g, K)
```

The functions `f₁` and `f₂` are highly nonlinear (10th order polynomials in the derivatives) - see `generate_equations.py` for the full symbolic forms.

### Boundary Conditions

**At anchors:** We specify either the anchor tension `T` or the natural length `n`. The angle is determined by shooting method (binary search on initial angle until the line reaches the target length/tension).

**At point masses:** Slackliners create discontinuities in the derivative fields. We use conservation of canonical momentum:

```
[∂L/∂y']_right - [∂L/∂y']_left = -M g    (weight causes momentum jump)
[∂L/∂n']_right - [∂L/∂n']_left = 0        (natural length continuous)
```

These jump conditions are solved numerically to propagate the solution across each point mass.

### Numerical Integration

The system is integrated using `scipy.integrate.RK45` with adaptive step sizing. For a given gap length and anchor tension:

1. **Standing tension solve:** Integrate the empty line (no slackliners) using a shooting method to find the initial anchor angle that produces the correct gap length at the target tension. This determines the natural length `n₀`.

2. **Loaded configuration:** If slackliners are present, re-integrate with natural length `n₀` fixed, applying jump boundary conditions at each mass location.

3. **Post-process:** Compute arc length `l`, tension `T = K(dl/dn - 1)`, and angle `A = arctan(|y'|)` from the solution.

The result is a physically accurate model of the slackline's 3D catenary shape under arbitrary loading conditions.

**Material properties** (Dyneemite Pro):
- `K = 40,000 N` (stiffness at 100% strain)
- `m = 0.115 kg/m` (linear mass density)
- `g = 9.81 m/s²`
