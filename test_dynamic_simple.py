#!/usr/bin/env python3
"""
Simple test of dynamic slackline simulation.
"""
import sys
sys.path.insert(0, '/home/user/slackline')

from src.api_dynamic import DynamicConstraints
import matplotlib
matplotlib.use('Agg')

print("Testing dynamic slackline simulation...")

# Create a simple 25m line with one person
print("\n1. Creating constraints...")
constraints = DynamicConstraints(gap_length=25, anchor_tension=2000,
                                n_nodes=50, damping_ratio=0.03)
constraints.add_slackliner(position=12.5, mass=75)
print("   ✓ Constraints created")

# Simulate a small pluck
print("\n2. Running simulation...")
dynamic_rig, static_rig = constraints.simulate_pluck(
    pluck_position=12.5,
    pluck_displacement=0.3,
    pluck_width=2.0,
    t_end=2.0,  # Short simulation
    n_frames=30  # Few frames for testing
)
print(f"   ✓ Simulation complete: {len(dynamic_rig.t)} frames")

# Save GIF
print("\n3. Generating GIF...")
dynamic_rig.save_gif('test_dynamic.gif', fps=15,
                    show_equilibrium=True,
                    equilibrium=static_rig.y)

print("\n✓ Test successful! Generated test_dynamic.gif")
