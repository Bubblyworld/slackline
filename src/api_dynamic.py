"""
High-level API for time-dependent slackline dynamics.

This module extends the static API to support dynamic simulations with GIF output.
"""
from src.api import Constraints, Rig
from src.core.dynamics import DynamicSlackline, create_slackliner_impulse, create_slackliner_oscillation, initial_pluck
from src.core.slacklines import dyneemite_pro
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from PIL import Image
import io


class DynamicRig:
    """
    A time-dependent rig that stores the evolution of the slackline.

    Attributes:
    -----------
    x : array
        Horizontal positions (fixed)
    t : array
        Time points
    y : array
        Vertical positions, shape (n_frames, n_nodes)
    """

    def __init__(self, x, t, y):
        """
        Initialize dynamic rig.

        Parameters:
        -----------
        x : array
            Horizontal positions
        t : array
            Time points
        y : array
            Vertical positions (n_frames, n_nodes)
        """
        self.x = x
        self.t = t
        self.y = y
        self.n_frames = len(t)
        self.n_nodes = len(x)

    def save_gif(self, filename, fps=30, dpi=100, figsize=(12, 6),
                 show_equilibrium=True, equilibrium=None):
        """
        Save animation as GIF.

        Parameters:
        -----------
        filename : str
            Output filename (should end in .gif)
        fps : int
            Frames per second
        dpi : int
            Resolution
        figsize : tuple
            Figure size (width, height) in inches
        show_equilibrium : bool
            Whether to show equilibrium line
        equilibrium : array or None
            Equilibrium y positions
        """
        print(f"Generating GIF with {self.n_frames} frames at {fps} fps...")

        fig, ax = plt.subplots(figsize=figsize)

        # Determine y-axis limits
        y_min = np.min(self.y) * 1.1
        y_max = max(0.1, np.max(self.y) * 1.1)

        # Plot equilibrium if provided
        if show_equilibrium and equilibrium is not None:
            ax.plot(self.x, equilibrium, 'k--', alpha=0.3, linewidth=1, label='Equilibrium')

        # Initialize line
        line, = ax.plot([], [], 'b-', linewidth=2, label='Slackline')
        time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                           fontsize=12, verticalalignment='top')

        ax.set_xlim(self.x[0], self.x[-1])
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel('Horizontal Distance (m)', fontsize=12)
        ax.set_ylabel('Vertical Position (m)', fontsize=12)
        ax.set_title('Dynamic Slackline Simulation', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')

        # Mark anchors
        ax.plot([self.x[0], self.x[-1]], [0, 0], 'ro', markersize=8, label='Anchors', zorder=5)

        def init():
            line.set_data([], [])
            time_text.set_text('')
            return line, time_text

        def animate(frame):
            line.set_data(self.x, self.y[frame])
            time_text.set_text(f't = {self.t[frame]:.3f} s')
            return line, time_text

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                      frames=self.n_frames, interval=1000/fps,
                                      blit=True, repeat=True)

        # Save as GIF
        print(f"Saving animation to {filename}...")
        writer = animation.PillowWriter(fps=fps)
        anim.save(filename, writer=writer, dpi=dpi)
        plt.close(fig)

        print(f"✓ GIF saved to {filename}")

    def save_plot_grid(self, filename, n_snapshots=6, figsize=(15, 10)):
        """
        Save a grid of snapshots at different times.

        Parameters:
        -----------
        filename : str
            Output filename
        n_snapshots : int
            Number of snapshots to show
        figsize : tuple
            Figure size
        """
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

        snapshot_indices = np.linspace(0, self.n_frames - 1, n_snapshots, dtype=int)

        y_min = np.min(self.y) * 1.1
        y_max = max(0.1, np.max(self.y) * 1.1)

        for idx, ax_idx in enumerate(snapshot_indices[:n_snapshots]):
            ax = axes[idx]
            ax.plot(self.x, self.y[ax_idx], 'b-', linewidth=2)
            ax.plot([self.x[0], self.x[-1]], [0, 0], 'ro', markersize=6)
            ax.set_xlim(self.x[0], self.x[-1])
            ax.set_ylim(y_min, y_max)
            ax.set_xlabel('Horizontal Distance (m)')
            ax.set_ylabel('Vertical Position (m)')
            ax.set_title(f't = {self.t[ax_idx]:.3f} s', fontweight='bold')
            ax.grid(True, alpha=0.3)

        plt.suptitle('Dynamic Slackline Snapshots', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"✓ Snapshot grid saved to {filename}")


class DynamicConstraints(Constraints):
    """
    Extension of Constraints for dynamic simulations.

    This class provides the same API as Constraints but adds methods for
    time-dependent simulation.
    """

    def __init__(self, slackline=dyneemite_pro, gap_length=100, anchor_tension=1000,
                 n_nodes=100, damping_ratio=0.02):
        """
        Initialize dynamic constraints.

        Parameters match Constraints, with additional:
        n_nodes : int
            Number of discretization nodes for dynamics
        damping_ratio : float
            Damping coefficient (fraction of critical damping)
        """
        super().__init__(slackline=slackline, gap_length=gap_length,
                        anchor_tension=anchor_tension)
        self.n_nodes = n_nodes
        self.damping_ratio = damping_ratio

    def simulate(self, t_span, initial_perturbation=None, slackliner_forces=None,
                n_frames=100):
        """
        Simulate time evolution from static equilibrium.

        Parameters:
        -----------
        t_span : tuple
            (t_start, t_end) time interval
        initial_perturbation : callable, array, or None
            Initial displacement from equilibrium
            - If callable: function f(x) -> displacement
            - If array: displacement values at each node
            - If None: start from equilibrium
        slackliner_forces : callable or None
            External forcing function f(t, x, y) -> force array
        n_frames : int
            Number of frames to output

        Returns:
        --------
        dynamic_rig : DynamicRig
            Time-dependent rig with x, t, y arrays
        """
        # First compute static equilibrium
        static_rig = self.rig()

        # Create dynamic system
        dyn = DynamicSlackline(static_rig, self.slackline,
                              n_nodes=self.n_nodes,
                              damping_ratio=self.damping_ratio)

        # Process initial perturbation
        if callable(initial_perturbation):
            pert = initial_perturbation(dyn.x)
        elif initial_perturbation is not None:
            pert = initial_perturbation
        else:
            pert = None

        # Run simulation
        t, y = dyn.simulate(t_span, slackliner_forces=slackliner_forces,
                           n_frames=n_frames, initial_perturbation=pert)

        # Create a simplified static rig for equilibrium reference (on same grid as dynamic)
        return DynamicRig(dyn.x, t, y), dyn.y_eq

    def simulate_pluck(self, pluck_position, pluck_displacement, t_end,
                      pluck_width=2.0, n_frames=100):
        """
        Simulate a pluck perturbation (pull and release).

        Parameters:
        -----------
        pluck_position : float
            Horizontal position of pluck
        pluck_displacement : float
            Upward displacement (meters)
        t_end : float
            Simulation duration (seconds)
        pluck_width : float
            Width of pluck (meters)
        n_frames : int
            Number of output frames

        Returns:
        --------
        dynamic_rig, static_rig : tuple
            Time-dependent and equilibrium rigs
        """
        perturbation = initial_pluck(pluck_position, pluck_displacement, pluck_width)
        return self.simulate((0, t_end), initial_perturbation=perturbation,
                           n_frames=n_frames)

    def simulate_bounce(self, bounce_position, frequency, amplitude, t_end,
                       phase=0, n_frames=100):
        """
        Simulate a slackliner bouncing at a fixed position.

        Parameters:
        -----------
        bounce_position : float
            Horizontal position of bouncing
        frequency : float
            Bounce frequency (Hz)
        amplitude : float
            Force amplitude (Newtons)
        t_end : float
            Simulation duration (seconds)
        phase : float
            Phase offset (radians)
        n_frames : int
            Number of output frames

        Returns:
        --------
        dynamic_rig, static_rig : tuple
            Time-dependent and equilibrium rigs
        """
        forces = create_slackliner_oscillation(bounce_position, frequency,
                                              amplitude, phase)
        return self.simulate((0, t_end), slackliner_forces=forces,
                           n_frames=n_frames)

    def simulate_impulse(self, impulse_position, impulse_magnitude,
                        impulse_duration, t_end, n_frames=100):
        """
        Simulate an impulse (jump) perturbation.

        Parameters:
        -----------
        impulse_position : float
            Horizontal position of impulse
        impulse_magnitude : float
            Force magnitude (Newtons)
        impulse_duration : float
            Duration of impulse (seconds)
        t_end : float
            Simulation duration (seconds)
        n_frames : int
            Number of output frames

        Returns:
        --------
        dynamic_rig, static_rig : tuple
            Time-dependent and equilibrium rigs
        """
        forces = create_slackliner_impulse(impulse_position, impulse_magnitude,
                                          impulse_duration)
        return self.simulate((0, t_end), slackliner_forces=forces,
                           n_frames=n_frames)
