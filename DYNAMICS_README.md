# Time-Dependent Slackline Dynamics

This extension implements **time-dependent** simulations of slackline physics using Lagrangian mechanics. While the original codebase solves for static equilibrium, this module computes the full dynamic evolution of the slackline, including wave propagation, oscillations, and response to time-varying forces.

## Overview

The dynamic solver discretizes the slackline into nodes and derives equations of motion from the Lagrangian including kinetic energy terms. The system evolves according to:

```
L = T - U = (kinetic energy) - (gravitational + elastic potential energy)
```

The equations of motion are integrated forward in time using `scipy.integrate.solve_ivp` with adaptive RK45 integration.

## Installation

The dynamic module uses the same dependencies as the static solver. No additional packages are required beyond the base installation:

```bash
poetry install
```

## Quick Start

### Basic Example: Pluck Simulation

```python
from src.api_dynamic import DynamicConstraints

# Create a 50m slackline with standing tension of 2500N
constraints = DynamicConstraints(gap_length=50, anchor_tension=2500)

# Add a 75kg person at the center
constraints.add_slackliner(position=25, mass=75)

# Simulate a pluck: pull center up 40cm and release
dynamic_rig, equilibrium = constraints.simulate_pluck(
    pluck_position=25.0,      # Position of perturbation
    pluck_displacement=0.4,    # 40cm upward
    pluck_width=2.0,          # 2m gaussian width
    t_end=5.0,                # 5 second simulation
    n_frames=100              # 100 frames output
)

# Save as animated GIF
dynamic_rig.save_gif('slackline_pluck.gif', fps=25,
                    show_equilibrium=True, equilibrium=equilibrium)
```

## API Reference

### `DynamicConstraints`

Main API class, extends `Constraints` from the static solver.

**Constructor:**
```python
DynamicConstraints(
    slackline=dyneemite_pro,  # Slackline parameters
    gap_length=100,            # Gap length (m)
    anchor_tension=1000,       # Standing tension (N)
    n_nodes=100,               # Number of discretization nodes
    damping_ratio=0.02         # Damping coefficient (fraction of critical)
)
```

**Parameters:**
- `slackline`: Slackline material properties (default: Dyneemite Pro)
- `gap_length`: Distance between anchors in meters
- `anchor_tension`: Standing anchor tension before loads applied
- `n_nodes`: Number of spatial discretization points (more = more accurate but slower)
- `damping_ratio`: Viscous damping coefficient (0-1, typical: 0.01-0.05)

### Simulation Methods

#### `simulate_pluck()`
Simulate a "pluck" perturbation - pull the line and release.

```python
dynamic_rig, equilibrium = constraints.simulate_pluck(
    pluck_position,      # Horizontal position of pluck (m)
    pluck_displacement,  # Vertical displacement (m, positive = up)
    t_end,              # Simulation duration (s)
    pluck_width=2.0,    # Gaussian width of perturbation (m)
    n_frames=100        # Number of output frames
)
```

**Returns:**
- `dynamic_rig`: DynamicRig object with time evolution
- `equilibrium`: Equilibrium configuration (array)

#### `simulate_bounce()`
Simulate a slackliner bouncing rhythmically.

```python
dynamic_rig, equilibrium = constraints.simulate_bounce(
    bounce_position,   # Horizontal position (m)
    frequency,         # Bounce frequency (Hz)
    amplitude,         # Force amplitude (N)
    t_end,            # Simulation duration (s)
    phase=0,          # Phase offset (radians)
    n_frames=100      # Number of output frames
)
```

#### `simulate_impulse()`
Simulate an impulsive force (e.g., a jump).

```python
dynamic_rig, equilibrium = constraints.simulate_impulse(
    impulse_position,   # Horizontal position (m)
    impulse_magnitude,  # Force magnitude (N, negative = downward)
    impulse_duration,   # Duration of impulse (s)
    t_end,             # Simulation duration (s)
    n_frames=100       # Number of output frames
)
```

#### `simulate()` (Advanced)
Low-level simulation with custom forcing.

```python
dynamic_rig, equilibrium = constraints.simulate(
    t_span,                 # (t_start, t_end)
    initial_perturbation,   # callable f(x) -> displacement or array
    slackliner_forces,      # callable f(t, x, y) -> force array
    n_frames=100           # Number of output frames
)
```

### `DynamicRig`

Output object containing time evolution data.

**Attributes:**
- `x`: Horizontal positions (m) - array shape (n_nodes,)
- `t`: Time points (s) - array shape (n_frames,)
- `y`: Vertical positions (m) - array shape (n_frames, n_nodes)

**Methods:**

```python
# Save as animated GIF
dynamic_rig.save_gif(
    filename,               # Output filename (e.g., 'anim.gif')
    fps=30,                # Frames per second
    dpi=100,               # Resolution
    figsize=(12, 6),       # Figure size (inches)
    show_equilibrium=True, # Show equilibrium line
    equilibrium=None       # Equilibrium y values
)

# Save snapshot grid
dynamic_rig.save_plot_grid(
    filename,           # Output filename (e.g., 'snapshots.png')
    n_snapshots=6,     # Number of snapshots to show
    figsize=(15, 10)   # Figure size (inches)
)
```

## Examples

### Example 1: Free Oscillation

```python
from src.api_dynamic import DynamicConstraints

# Empty 50m line
constraints = DynamicConstraints(gap_length=50, anchor_tension=2000, n_nodes=100)

# Pluck asymmetrically
dynamic_rig, eq = constraints.simulate_pluck(
    pluck_position=20.0,    # Off-center
    pluck_displacement=0.4,
    t_end=8.0,
    n_frames=200
)

dynamic_rig.save_gif('free_oscillation.gif', fps=30, equilibrium=eq)
```

### Example 2: Person Bouncing

```python
from src.api_dynamic import DynamicConstraints

# 100m line with 80kg person
constraints = DynamicConstraints(gap_length=100, anchor_tension=3000, n_nodes=120)
constraints.add_slackliner(position=40, mass=80)

# Person bounces at 1 Hz
dynamic_rig, eq = constraints.simulate_bounce(
    bounce_position=40.0,
    frequency=1.0,         # 1 Hz
    amplitude=500.0,       # 500N force amplitude
    t_end=8.0,
    n_frames=200
)

dynamic_rig.save_gif('person_bouncing.gif', fps=30, equilibrium=eq)
```

### Example 3: Jump Landing

```python
from src.api_dynamic import DynamicConstraints

# 25m line with 70kg person
constraints = DynamicConstraints(gap_length=25, anchor_tension=2000, n_nodes=80)
constraints.add_slackliner(position=12.5, mass=70)

# Sudden downward impulse (landing)
dynamic_rig, eq = constraints.simulate_impulse(
    impulse_position=12.5,
    impulse_magnitude=-800.0,  # 800N downward
    impulse_duration=0.2,      # 0.2s impact
    t_end=6.0,
    n_frames=150
)

dynamic_rig.save_gif('jump_landing.gif', fps=30, equilibrium=eq)
```

### Example 4: Two People

```python
from src.api_dynamic import DynamicConstraints

# 100m line with two people
constraints = DynamicConstraints(gap_length=100, anchor_tension=3500, n_nodes=150)
constraints.add_slackliner(position=30, mass=70)
constraints.add_slackliner(position=70, mass=80)

# Pluck between them
dynamic_rig, eq = constraints.simulate_pluck(
    pluck_position=50.0,
    pluck_displacement=0.3,
    t_end=12.0,
    n_frames=240
)

dynamic_rig.save_gif('two_people.gif', fps=30, equilibrium=eq)
```

## Physics Details

### Discretization

The continuous slackline is discretized into N nodes with positions `(x_i, y_i(t))` where:
- `x_i` are fixed horizontal positions
- `y_i(t)` are time-dependent vertical positions

### Equations of Motion

For each internal node i, the equation of motion is:

```
m_i * ÿ_i = F_elastic,i + F_gravity,i + F_damping,i + F_external,i
```

Where:
- **Elastic forces** come from tension in adjacent segments
- **Gravity** acts downward: `F_g = -m_i * g`
- **Damping** is viscous: `F_d = -c * ẏ_i`
- **External forces** include slackliner inputs

### Elastic Tension

Each segment [i, i+1] has:
- Current length: `l = √(Δx² + Δy²)`
- Natural length: `n` (from equilibrium)
- Strain: `ε = (l - n) / n`
- Tension: `T = K * ε`

The force on node i from this segment is `F = T * sin(θ)` where `θ` is the angle.

### Damping

Viscous damping with coefficient:
```
c = ζ * 2 * √(K * m / n)
```
where `ζ` is the damping ratio (typically 0.01-0.05).

### Numerical Integration

The ODE system is integrated using `scipy.integrate.solve_ivp` with:
- Method: RK45 (adaptive 4th/5th order Runge-Kutta)
- Relative tolerance: 1e-6
- Absolute tolerance: 1e-8

## Performance Tips

1. **Fewer nodes = faster**: Start with `n_nodes=50-80` for testing
2. **Shorter simulations**: Use `t_end=2-5` seconds initially
3. **Fewer frames**: Use `n_frames=50-100` for quick previews
4. **Damping**: Higher damping (0.03-0.05) stabilizes faster

## Limitations

1. **Small deflections**: Works best for deflections < ~20% of gap length
2. **No buckling**: Compression forces are limited to prevent numerical instability
3. **Fixed horizontal positions**: Nodes don't move horizontally (valid for small deflections)
4. **Point masses**: Slackliners are treated as point masses, not distributed loads

## Advanced: Custom Forcing

For custom time-dependent forces, use the `simulate()` method:

```python
def my_force_function(t, x, y):
    """
    Custom forcing function.

    Args:
        t: current time (float)
        x: node positions (array)
        y: current vertical positions (array)

    Returns:
        forces: force array (same shape as x)
    """
    forces = np.zeros_like(x)

    # Apply force at x=25m with time variation
    idx = np.argmin(np.abs(x - 25.0))
    forces[idx] = 100 * np.sin(2 * np.pi * t)  # Sinusoidal force

    return forces

# Use custom forcing
dynamic_rig, eq = constraints.simulate(
    t_span=(0, 10),
    slackliner_forces=my_force_function,
    n_frames=200
)
```

## Running the Examples

A comprehensive example script is provided:

```bash
poetry run python generate_dynamic_examples.py
```

This generates 5 different scenarios with GIF animations and snapshot grids.

## Theory: From Static to Dynamic

### Static Lagrangian

The static solver uses:
```
L = m*g*y*n' + K/2*(1 + y'²)/n' - K*√(1 + y'²) + K/2*n'
```

### Dynamic Extension

For dynamics, we add kinetic energy:
```
T = Σ (1/2) * m_i * (ẏ_i)²
```

The full Lagrangian becomes:
```
L = T - U = (kinetic energy) - (potential energy)
```

Applying the Euler-Lagrange equations:
```
d/dt(∂L/∂ẏ_i) - ∂L/∂y_i = 0
```

yields the equations of motion for each node.

## Troubleshooting

**"Simulation complete. Generated 1 frames"**
- Solution: Increase integration time or check for numerical instabilities

**Warnings about overflow**
- Solution: Increase `damping_ratio` or decrease `n_nodes`

**GIF looks choppy**
- Solution: Increase `n_frames` (100-300 recommended) and `fps` (25-60)

**Simulation too slow**
- Solution: Decrease `n_nodes` (50-100 is usually sufficient) or `n_frames`

## Citation

If you use this dynamic extension in your research, please cite the original static solver and mention the time-dependent extension.

## License

MIT License (same as parent project)
