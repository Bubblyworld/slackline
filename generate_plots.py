#!/usr/bin/env python3
"""
Script to generate example slackline plots for different configurations.
"""
import sys
sys.path.insert(0, '/home/user/slackline')

from src.api import Constraints, Rig
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

def generate_plot_1():
    """25m slackline with 80kg person in the middle"""
    print("Generating Plot 1: 25m slackline with 80kg person at center...")
    constraints = Constraints(gap_length=25, anchor_tension=2000)
    constraints.add_slackliner(position=12.5, mass=80)
    rig = constraints.rig()

    plt.figure(figsize=(12, 10))

    plt.subplot(221)
    plt.plot(rig.x, rig.y, 'b-', linewidth=2)
    plt.axvline(x=12.5, color='r', linestyle='--', alpha=0.7, label='80kg')
    plt.title("Curve - 25m Line", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Vertical Drop (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(222)
    plt.plot(rig.x, rig.n, label="Natural Length", linewidth=2)
    plt.plot(rig.x, rig.l, label="Arc Length", linewidth=2)
    plt.axvline(x=12.5, color='r', linestyle='--', alpha=0.5)
    plt.title("Natural and Arc Lengths", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Length (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(223)
    plt.plot(rig.x, rig.T, 'r-', linewidth=2)
    plt.axvline(x=12.5, color='r', linestyle='--', alpha=0.5)
    plt.title("Tension Distribution", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Tension (N)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(224)
    plt.plot(rig.x, rig.A, 'g-', linewidth=2)
    plt.axvline(x=12.5, color='r', linestyle='--', alpha=0.5)
    plt.title("Angle from Horizontal", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Angle (degrees)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/slackline/plot1_25m.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved as plot1_25m.png")

def generate_plot_2():
    """100m slackline with two people"""
    print("Generating Plot 2: 100m slackline with two people...")
    constraints = Constraints(gap_length=100, anchor_tension=3000)
    constraints.add_slackliner(position=30, mass=70)
    constraints.add_slackliner(position=70, mass=80)
    rig = constraints.rig()

    plt.figure(figsize=(12, 10))

    plt.subplot(221)
    plt.plot(rig.x, rig.y, 'b-', linewidth=2)
    plt.axvline(x=30, color='r', linestyle='--', alpha=0.7, label='70kg')
    plt.axvline(x=70, color='orange', linestyle='--', alpha=0.7, label='80kg')
    plt.title("Curve - 100m Line", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Vertical Drop (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(222)
    plt.plot(rig.x, rig.n, label="Natural Length", linewidth=2)
    plt.plot(rig.x, rig.l, label="Arc Length", linewidth=2)
    plt.axvline(x=30, color='r', linestyle='--', alpha=0.5)
    plt.axvline(x=70, color='orange', linestyle='--', alpha=0.5)
    plt.title("Natural and Arc Lengths", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Length (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(223)
    plt.plot(rig.x, rig.T, 'r-', linewidth=2)
    plt.axvline(x=30, color='r', linestyle='--', alpha=0.5)
    plt.axvline(x=70, color='orange', linestyle='--', alpha=0.5)
    plt.title("Tension Distribution", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Tension (N)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(224)
    plt.plot(rig.x, rig.A, 'g-', linewidth=2)
    plt.axvline(x=30, color='r', linestyle='--', alpha=0.5)
    plt.axvline(x=70, color='orange', linestyle='--', alpha=0.5)
    plt.title("Angle from Horizontal", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Angle (degrees)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/slackline/plot2_100m.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved as plot2_100m.png")

def generate_plot_3():
    """500m slackline with someone 100m out"""
    print("Generating Plot 3: 500m slackline with person at 100m...")
    constraints = Constraints(gap_length=500, anchor_tension=8000)
    constraints.add_slackliner(position=100, mass=75)
    rig = constraints.rig()

    plt.figure(figsize=(12, 10))

    plt.subplot(221)
    plt.plot(rig.x, rig.y, 'b-', linewidth=2)
    plt.axvline(x=100, color='r', linestyle='--', alpha=0.7, label='75kg')
    plt.title("Curve - 500m Line", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Vertical Drop (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(222)
    plt.plot(rig.x, rig.n, label="Natural Length", linewidth=2)
    plt.plot(rig.x, rig.l, label="Arc Length", linewidth=2)
    plt.axvline(x=100, color='r', linestyle='--', alpha=0.5)
    plt.title("Natural and Arc Lengths", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Length (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(223)
    plt.plot(rig.x, rig.T, 'r-', linewidth=2)
    plt.axvline(x=100, color='r', linestyle='--', alpha=0.5)
    plt.title("Tension Distribution", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Tension (N)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(224)
    plt.plot(rig.x, rig.A, 'g-', linewidth=2)
    plt.axvline(x=100, color='r', linestyle='--', alpha=0.5)
    plt.title("Angle from Horizontal", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Angle (degrees)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/slackline/plot3_500m.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved as plot3_500m.png")


if __name__ == "__main__":
    print("\n=== Generating Slackline Plots ===\n")
    generate_plot_1()
    generate_plot_2()
    generate_plot_3()
    print("\n=== All plots generated successfully! ===\n")
