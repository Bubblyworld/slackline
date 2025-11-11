#!/usr/bin/env python3
"""
Generate example dynamic slackline simulations with GIF outputs.

This script demonstrates the time-dependent Euler-Lagrange solver for
slackline dynamics, showing various scenarios like plucks, bounces, and impulses.
"""
import sys
sys.path.insert(0, '/home/user/slackline')

from src.api_dynamic import DynamicConstraints
from src.core.dynamics import initial_pluck
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend


def example_1_pluck():
    """
    Example 1: 50m slackline with 75kg person, plucked upward at center.

    This simulates pulling the center of the line upward and releasing,
    showing the resulting wave propagation and oscillations.
    """
    print("\n" + "="*70)
    print("Example 1: Pluck Perturbation (50m line, 75kg person at center)")
    print("="*70)

    # Create constraints with a person at the center
    constraints = DynamicConstraints(gap_length=50, anchor_tension=2500,
                                    n_nodes=100, damping_ratio=0.03)
    constraints.add_slackliner(position=25, mass=75)

    # Simulate a pluck: pull center up 0.5m and release
    dynamic_rig, static_rig = constraints.simulate_pluck(
        pluck_position=25.0,    # At the center
        pluck_displacement=0.5,  # 50cm upward
        pluck_width=2.0,        # 2m wide gaussian
        t_end=10.0,             # 10 second simulation
        n_frames=200            # 200 frames = smooth animation
    )

    # Save outputs
    dynamic_rig.save_gif('dynamic_pluck_50m.gif', fps=30,
                        show_equilibrium=True,
                        equilibrium=static_rig.y)
    dynamic_rig.save_plot_grid('dynamic_pluck_50m_snapshots.png', n_snapshots=6)

    print(f"✓ Example 1 complete\n")


def example_2_bounce():
    """
    Example 2: 100m slackline with person bouncing at 1 Hz.

    This simulates a slackliner bouncing rhythmically, creating periodic
    oscillations that propagate along the line.
    """
    print("\n" + "="*70)
    print("Example 2: Bouncing Slackliner (100m line, 1 Hz bounce)")
    print("="*70)

    # Create constraints
    constraints = DynamicConstraints(gap_length=100, anchor_tension=3000,
                                    n_nodes=120, damping_ratio=0.02)
    constraints.add_slackliner(position=40, mass=80)

    # Simulate bouncing: 1 Hz frequency, 500N amplitude
    dynamic_rig, static_rig = constraints.simulate_bounce(
        bounce_position=40.0,   # 40m from left anchor
        frequency=1.0,          # 1 Hz bouncing
        amplitude=500.0,        # 500N force amplitude
        t_end=8.0,              # 8 second simulation
        n_frames=200            # 200 frames
    )

    # Save outputs
    dynamic_rig.save_gif('dynamic_bounce_100m.gif', fps=30,
                        show_equilibrium=True,
                        equilibrium=static_rig.y)
    dynamic_rig.save_plot_grid('dynamic_bounce_100m_snapshots.png', n_snapshots=6)

    print(f"✓ Example 2 complete\n")


def example_3_impulse():
    """
    Example 3: 25m slackline with impulsive jump.

    This simulates a sudden jump (impulse force), showing the shock wave
    propagation and subsequent oscillations.
    """
    print("\n" + "="*70)
    print("Example 3: Impulse/Jump (25m line, 70kg person)")
    print("="*70)

    # Create constraints
    constraints = DynamicConstraints(gap_length=25, anchor_tension=2000,
                                    n_nodes=80, damping_ratio=0.025)
    constraints.add_slackliner(position=12.5, mass=70)

    # Simulate an impulse: strong downward force for 0.2 seconds
    dynamic_rig, static_rig = constraints.simulate_impulse(
        impulse_position=12.5,   # At center
        impulse_magnitude=-800.0, # 800N downward (negative = down)
        impulse_duration=0.2,    # 0.2 second impulse
        t_end=6.0,               # 6 second simulation
        n_frames=150             # 150 frames
    )

    # Save outputs
    dynamic_rig.save_gif('dynamic_impulse_25m.gif', fps=30,
                        show_equilibrium=True,
                        equilibrium=static_rig.y)
    dynamic_rig.save_plot_grid('dynamic_impulse_25m_snapshots.png', n_snapshots=6)

    print(f"✓ Example 3 complete\n")


def example_4_two_people():
    """
    Example 4: 100m line with two people, complex initial perturbation.

    This shows interactions between two masses and more complex wave dynamics.
    """
    print("\n" + "="*70)
    print("Example 4: Two People with Pluck (100m line)")
    print("="*70)

    # Create constraints with two people
    constraints = DynamicConstraints(gap_length=100, anchor_tension=3500,
                                    n_nodes=150, damping_ratio=0.02)
    constraints.add_slackliner(position=30, mass=70)
    constraints.add_slackliner(position=70, mass=80)

    # Simulate a pluck between the two people
    dynamic_rig, static_rig = constraints.simulate_pluck(
        pluck_position=50.0,     # Midpoint between them
        pluck_displacement=0.3,   # 30cm upward
        pluck_width=3.0,         # 3m wide
        t_end=12.0,              # 12 second simulation
        n_frames=240             # 240 frames for smooth motion
    )

    # Save outputs
    dynamic_rig.save_gif('dynamic_two_people_100m.gif', fps=30,
                        show_equilibrium=True,
                        equilibrium=static_rig.y)
    dynamic_rig.save_plot_grid('dynamic_two_people_100m_snapshots.png', n_snapshots=6)

    print(f"✓ Example 4 complete\n")


def example_5_free_oscillation():
    """
    Example 5: Free oscillation of empty line.

    This shows the natural modes of vibration of an unloaded slackline.
    """
    print("\n" + "="*70)
    print("Example 5: Free Oscillation (50m empty line)")
    print("="*70)

    # Create constraints with NO slackliners
    constraints = DynamicConstraints(gap_length=50, anchor_tension=2000,
                                    n_nodes=100, damping_ratio=0.01)

    # Simulate with a simple pluck
    dynamic_rig, static_rig = constraints.simulate_pluck(
        pluck_position=20.0,     # Asymmetric position
        pluck_displacement=0.4,   # 40cm upward
        pluck_width=2.5,         # 2.5m wide
        t_end=8.0,               # 8 second simulation
        n_frames=200             # 200 frames
    )

    # Save outputs
    dynamic_rig.save_gif('dynamic_free_oscillation_50m.gif', fps=30,
                        show_equilibrium=True,
                        equilibrium=static_rig.y)
    dynamic_rig.save_plot_grid('dynamic_free_oscillation_50m_snapshots.png',
                              n_snapshots=6)

    print(f"✓ Example 5 complete\n")


if __name__ == "__main__":
    print("\n" + "="*70)
    print("DYNAMIC SLACKLINE SIMULATION EXAMPLES")
    print("Time-dependent Euler-Lagrange Equations")
    print("="*70)

    # Run all examples
    try:
        example_1_pluck()
        example_2_bounce()
        example_3_impulse()
        example_4_two_people()
        example_5_free_oscillation()

        print("\n" + "="*70)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("="*70)
        print("\nGenerated files:")
        print("  - dynamic_pluck_50m.gif (+ snapshots)")
        print("  - dynamic_bounce_100m.gif (+ snapshots)")
        print("  - dynamic_impulse_25m.gif (+ snapshots)")
        print("  - dynamic_two_people_100m.gif (+ snapshots)")
        print("  - dynamic_free_oscillation_50m.gif (+ snapshots)")
        print("\n")

    except Exception as e:
        print(f"\n❌ Error during simulation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
