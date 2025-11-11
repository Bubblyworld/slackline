# Dynamic Slackline Simulations - Quick Start

## What's New: Time-Dependent Dynamics

This extension adds **time-dependent** simulations to the static slackline solver. You can now:

✨ Simulate wave propagation and oscillations
✨ Model slackliners bouncing and jumping
✨ Generate animated GIFs of dynamic behavior
✨ Study transient responses to perturbations

## 30-Second Example

```python
from src.api_dynamic import DynamicConstraints

# Create 50m slackline with person
constraints = DynamicConstraints(gap_length=50, anchor_tension=2500)
constraints.add_slackliner(position=25, mass=75)

# Simulate pluck and save as GIF
dynamic_rig, equilibrium = constraints.simulate_pluck(
    pluck_position=25, pluck_displacement=0.4, t_end=5.0, n_frames=100
)
dynamic_rig.save_gif('output.gif', fps=25, equilibrium=equilibrium)
```

That's it! You'll get an animated GIF showing the slackline oscillating over 5 seconds.

## Key Features

### Same API as Static Solver
```python
# Static equilibrium (existing)
from src.api import Constraints
constraints = Constraints(gap_length=50, anchor_tension=2000)
constraints.add_slackliner(position=25, mass=75)
rig = constraints.rig()  # Returns static solution

# Dynamic evolution (new!)
from src.api_dynamic import DynamicConstraints
constraints = DynamicConstraints(gap_length=50, anchor_tension=2000)
constraints.add_slackliner(position=25, mass=75)
dynamic_rig, eq = constraints.simulate_pluck(...)  # Returns time series
```

### Three Built-in Scenarios

1. **Pluck**: Pull the line and release
```python
dynamic_rig, eq = constraints.simulate_pluck(
    pluck_position=25, pluck_displacement=0.4, t_end=5.0
)
```

2. **Bounce**: Rhythmic bouncing
```python
dynamic_rig, eq = constraints.simulate_bounce(
    bounce_position=25, frequency=1.0, amplitude=500, t_end=8.0
)
```

3. **Impulse**: Jump or sudden force
```python
dynamic_rig, eq = constraints.simulate_impulse(
    impulse_position=25, impulse_magnitude=-800, impulse_duration=0.2, t_end=6.0
)
```

### GIF Output

```python
dynamic_rig.save_gif(
    'animation.gif',
    fps=30,                 # Frames per second
    show_equilibrium=True,  # Show static equilibrium line
    equilibrium=eq          # Equilibrium configuration
)
```

## How It Works

### Physics
- Extends Lagrangian mechanics to include **kinetic energy**
- Discretizes slackline into nodes with equations of motion:
  ```
  m_i * ÿ_i = F_elastic + F_gravity + F_damping + F_external
  ```
- Integrates forward in time using adaptive RK45 method

### Parameters
```python
DynamicConstraints(
    gap_length=100,      # Same as static solver
    anchor_tension=1000, # Same as static solver
    n_nodes=100,         # Discretization resolution (50-150 typical)
    damping_ratio=0.02   # Viscous damping (0.01-0.05 typical)
)
```

## Running Examples

```bash
# Run comprehensive examples
poetry run python generate_dynamic_examples.py

# Run quick test
poetry run python test_dynamic_minimal.py
```

This generates multiple GIF animations showing different dynamic scenarios.

## Documentation

See **`DYNAMICS_README.md`** for:
- Complete API reference
- Physics details and equations
- Performance tips
- Troubleshooting guide
- Advanced custom forcing examples

## Files Added

```
src/core/dynamics.py           # Core dynamics implementation
src/api_dynamic.py             # High-level API
generate_dynamic_examples.py   # Example scenarios
test_dynamic_minimal.py        # Simple test
test_dynamic_simple.py         # Test with person
DYNAMICS_README.md             # Full documentation
DYNAMIC_QUICKSTART.md          # This file
```

## Performance

Typical simulation (50m, 100 nodes, 5 seconds, 100 frames):
- Computation time: ~10-30 seconds
- Output GIF size: ~2-5 MB

Tips for faster simulations:
- Fewer nodes: `n_nodes=50-80`
- Shorter time: `t_end=2-3`
- Fewer frames: `n_frames=50`

## Requirements

No additional dependencies beyond the static solver:
- scipy (ODE integration)
- numpy (arrays)
- matplotlib (plotting/animation)
- pillow (GIF generation)

All already included in `pyproject.toml`.
