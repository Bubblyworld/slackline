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
    """Short slackline (25m) with no slackliners - beginner setup"""
    print("Generating Plot 1: Short slackline (25m, 1500N)...")
    constraints = Constraints(gap_length=25, anchor_tension=1500)
    rig = constraints.rig()

    plt.figure(figsize=(12, 10))

    plt.subplot(221)
    plt.plot(rig.x, rig.y, 'b-', linewidth=2)
    plt.title("Curve - 25m Beginner Line", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Vertical Drop (m)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(222)
    plt.plot(rig.x, rig.n, label="Natural Length", linewidth=2)
    plt.plot(rig.x, rig.l, label="Arc Length", linewidth=2)
    plt.title("Natural and Arc Lengths", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Length (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(223)
    plt.plot(rig.x, rig.T, 'r-', linewidth=2)
    plt.title("Tension Distribution", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Tension (N)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(224)
    plt.plot(rig.x, rig.A, 'g-', linewidth=2)
    plt.title("Angle from Horizontal", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Angle (degrees)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/slackline/plot1_short_line.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved as plot1_short_line.png")

def generate_plot_2():
    """Medium slackline (50m) with one person in the middle"""
    print("Generating Plot 2: Medium slackline (50m) with 75kg slackliner at center...")
    constraints = Constraints(gap_length=50, anchor_tension=2000)
    constraints.add_slackliner(position=25, mass=75)
    rig = constraints.rig()

    plt.figure(figsize=(12, 10))

    plt.subplot(221)
    plt.plot(rig.x, rig.y, 'b-', linewidth=2)
    plt.axvline(x=25, color='r', linestyle='--', alpha=0.7, label='Slackliner (75kg)')
    plt.title("Curve - 50m Line with Person at Center", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Vertical Drop (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(222)
    plt.plot(rig.x, rig.n, label="Natural Length", linewidth=2)
    plt.plot(rig.x, rig.l, label="Arc Length", linewidth=2)
    plt.axvline(x=25, color='r', linestyle='--', alpha=0.5)
    plt.title("Natural and Arc Lengths", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Length (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(223)
    plt.plot(rig.x, rig.T, 'r-', linewidth=2)
    plt.axvline(x=25, color='r', linestyle='--', alpha=0.5)
    plt.title("Tension Distribution", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Tension (N)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(224)
    plt.plot(rig.x, rig.A, 'g-', linewidth=2)
    plt.axvline(x=25, color='r', linestyle='--', alpha=0.5)
    plt.title("Angle from Horizontal", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Angle (degrees)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/slackline/plot2_medium_line_one_person.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved as plot2_medium_line_one_person.png")

def generate_plot_3():
    """Long slackline (100m) with two people"""
    print("Generating Plot 3: Long slackline (100m) with two slackliners...")
    constraints = Constraints(gap_length=100, anchor_tension=3000)
    constraints.add_slackliner(position=30, mass=70)
    constraints.add_slackliner(position=70, mass=80)
    rig = constraints.rig()

    plt.figure(figsize=(12, 10))

    plt.subplot(221)
    plt.plot(rig.x, rig.y, 'b-', linewidth=2)
    plt.axvline(x=30, color='r', linestyle='--', alpha=0.7, label='Person 1 (70kg)')
    plt.axvline(x=70, color='orange', linestyle='--', alpha=0.7, label='Person 2 (80kg)')
    plt.title("Curve - 100m Line with Two People", fontsize=14, fontweight='bold')
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
    plt.savefig('/home/user/slackline/plot3_long_line_two_people.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved as plot3_long_line_two_people.png")

def generate_plot_4():
    """Very long highline (200m) - extreme case"""
    print("Generating Plot 4: Very long highline (200m, high tension)...")
    constraints = Constraints(gap_length=200, anchor_tension=5000)
    rig = constraints.rig()

    plt.figure(figsize=(12, 10))

    plt.subplot(221)
    plt.plot(rig.x, rig.y, 'b-', linewidth=2)
    plt.title("Curve - 200m Highline (Empty)", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Vertical Drop (m)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(222)
    plt.plot(rig.x, rig.n, label="Natural Length", linewidth=2)
    plt.plot(rig.x, rig.l, label="Arc Length", linewidth=2)
    plt.title("Natural and Arc Lengths", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Length (m)", fontsize=11)
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(223)
    plt.plot(rig.x, rig.T, 'r-', linewidth=2)
    plt.title("Tension Distribution", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Tension (N)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.subplot(224)
    plt.plot(rig.x, rig.A, 'g-', linewidth=2)
    plt.title("Angle from Horizontal", fontsize=14, fontweight='bold')
    plt.xlabel("Horizontal Distance (m)", fontsize=11)
    plt.ylabel("Angle (degrees)", fontsize=11)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/slackline/plot4_extreme_highline.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved as plot4_extreme_highline.png")

if __name__ == "__main__":
    print("\n=== Generating Slackline Plots ===\n")
    generate_plot_1()
    generate_plot_2()
    generate_plot_3()
    generate_plot_4()
    print("\n=== All plots generated successfully! ===\n")
