#!/usr/bin/env python3
"""
Minimal test of dynamic slackline simulation - no slackliners, just pluck.
"""
import sys
sys.path.insert(0, '/home/user/slackline')

from src.api_dynamic import DynamicConstraints
import matplotlib
matplotlib.use('Agg')

print("Testing minimal dynamic slackline simulation...")

# Create a simple 25m line with NO slackliners
print("\n1. Creating constraints...")
constraints = DynamicConstraints(gap_length=25, anchor_tension=2000,
                                n_nodes=30, damping_ratio=0.05)
print("   ✓ Constraints created")

# Simulate a small pluck
print("\n2. Running simulation...")
dynamic_rig, equilibrium = constraints.simulate_pluck(
    pluck_position=12.5,
    pluck_displacement=0.2,
    pluck_width=1.5,
    t_end=1.5,  # Very short simulation
    n_frames=20  # Very few frames
)
print(f"   ✓ Simulation complete: {len(dynamic_rig.t)} frames")

# Save GIF
print("\n3. Generating GIF...")
dynamic_rig.save_gif('test_minimal.gif', fps=10,
                    show_equilibrium=True,
                    equilibrium=equilibrium)

print("\n✓ Minimal test successful! Generated test_minimal.gif")
