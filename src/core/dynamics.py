"""
Time-dependent slackline dynamics using Lagrangian mechanics.

This module extends the static equilibrium solver to handle dynamic evolution
of the slackline under time-varying loads. The approach discretizes the slackline
into nodes and derives equations of motion from the Lagrangian including kinetic
energy terms.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d


class DynamicSlackline:
    """
    Represents a discretized slackline for dynamic time evolution.

    The slackline is divided into N nodes, each with position (x_i, y_i(t)).
    The horizontal positions x_i are fixed (small horizontal motion neglected),
    while vertical positions y_i(t) evolve according to the equations of motion.
    """

    def __init__(self, rig, slackline, n_nodes=100, damping_ratio=0.02):
        """
        Initialize dynamic slackline from static equilibrium.

        Parameters:
        -----------
        rig : Rig
            Static equilibrium configuration
        slackline : Slackline
            Slackline parameters (m, g, K)
        n_nodes : int
            Number of discretization nodes (default: 100)
        damping_ratio : float
            Damping coefficient (fraction of critical damping, default: 0.02)
        """
        self.slackline = slackline
        self.n_nodes = n_nodes
        self.damping_ratio = damping_ratio

        # Discretize the equilibrium configuration
        self.x = np.linspace(rig.x[0], rig.x[-1], n_nodes)

        # Interpolate static solution to node positions
        self.y_eq = np.interp(self.x, rig.x, rig.y)
        self.n_eq = np.interp(self.x, rig.x, rig.n)

        # Compute segment properties
        self.dx = self.x[1] - self.x[0]

        # Store equilibrium arc length and derivatives for each segment
        y_x_eq = np.gradient(self.y_eq, self.x)
        self.dl_eq = np.sqrt(1 + y_x_eq**2) * self.dx

        # Segment masses (mass per unit arc length * arc length)
        self.segment_mass = slackline.m * self.dl_eq

        # Natural lengths from equilibrium
        self.dn_eq = np.gradient(self.n_eq, self.x) * self.dx

    def kinetic_energy(self, y_dot):
        """
        Compute kinetic energy of the system.

        T = Σ (1/2) * m_i * v_i²
        """
        return 0.5 * np.sum(self.segment_mass * y_dot**2)

    def potential_energy(self, y):
        """
        Compute potential energy (gravitational + elastic).

        U = U_gravity + U_elastic
        """
        # Gravitational potential energy
        U_grav = np.sum(self.segment_mass * self.slackline.g * y)

        # Elastic potential energy
        # Compute current arc length
        y_x = np.gradient(y, self.x)
        dl = np.sqrt(1 + y_x**2) * self.dx

        # Strain in each segment
        epsilon = (dl - self.dn_eq) / self.dn_eq

        # Elastic energy: (1/2) K ε² per segment
        U_elastic = 0.5 * self.slackline.K * np.sum(epsilon**2 * self.dn_eq)

        return U_grav + U_elastic

    def equations_of_motion(self, t, state, slackliner_forces=None):
        """
        Compute time derivatives for the system state.

        State vector: [y_1, ..., y_N, ẏ_1, ..., ẏ_N]

        Returns: [ẏ_1, ..., ẏ_N, ÿ_1, ..., ÿ_N]

        Parameters:
        -----------
        t : float
            Current time
        state : array
            Current state [positions, velocities]
        slackliner_forces : callable or None
            Function f(t, x, y) -> force array for external forcing
        """
        n = self.n_nodes
        y = state[:n]
        y_dot = state[n:]

        # Compute forces from elastic potential energy
        # F_i = -∂U/∂y_i

        # Compute current configuration
        y_x = np.gradient(y, self.x)
        dl = np.sqrt(1 + y_x**2) * self.dx

        # Tension in each segment: T = K * (dl/dn - 1)
        tension = self.slackline.K * (dl / self.dn_eq - 1)

        # Force from tension (using finite differences)
        # F_y,i = -∂/∂y_i [ Σ (1/2) T_j (dl_j)² ]
        # This requires computing how each node affects adjacent segments

        forces = np.zeros(n)

        for i in range(n):
            if i == 0 or i == n-1:
                # Boundary nodes are fixed
                forces[i] = 0
                continue

            # Node i affects segments [i-1, i] and [i, i+1]
            # Compute force from elastic tension

            # Segment to the left [i-1, i]
            dy_left = y[i] - y[i-1]
            dx_left = self.dx
            dl_left = np.sqrt(dx_left**2 + dy_left**2)
            # Add numerical safeguards
            if dl_left > 1e-10 and np.isfinite(dl_left):
                strain_left = (dl_left - self.dn_eq[i-1]) / self.dn_eq[i-1]
                # Limit extreme strains to prevent numerical overflow
                strain_left = np.clip(strain_left, -0.5, 2.0)
                T_left = self.slackline.K * strain_left
                # Force component in y direction
                sin_theta_left = dy_left / dl_left
                forces[i] += T_left * sin_theta_left

            # Segment to the right [i, i+1]
            dy_right = y[i+1] - y[i]
            dx_right = self.dx
            dl_right = np.sqrt(dx_right**2 + dy_right**2)
            # Add numerical safeguards
            if dl_right > 1e-10 and np.isfinite(dl_right):
                strain_right = (dl_right - self.dn_eq[i]) / self.dn_eq[i]
                # Limit extreme strains to prevent numerical overflow
                strain_right = np.clip(strain_right, -0.5, 2.0)
                T_right = self.slackline.K * strain_right
                # Force component in y direction (pointing right, so negative contribution)
                sin_theta_right = dy_right / dl_right
                forces[i] -= T_right * sin_theta_right

        # Gravitational force
        forces -= self.segment_mass * self.slackline.g

        # Damping force: F_damp = -c * v
        # Use viscous damping proportional to velocity
        damping_coeff = self.damping_ratio * 2 * np.sqrt(self.slackline.K * self.segment_mass / self.dn_eq)
        forces -= damping_coeff * y_dot

        # External forces (e.g., from slackliners jumping)
        if slackliner_forces is not None:
            forces += slackliner_forces(t, self.x, y)

        # Accelerations: a = F/m
        y_ddot = forces / self.segment_mass

        # Boundary conditions: anchors are fixed
        y_ddot[0] = 0
        y_ddot[-1] = 0

        # Return derivative of state
        return np.concatenate([y_dot, y_ddot])

    def simulate(self, t_span, slackliner_forces=None, n_frames=100,
                 initial_perturbation=None):
        """
        Simulate time evolution of the slackline.

        Parameters:
        -----------
        t_span : tuple
            (t_start, t_end) time interval
        slackliner_forces : callable or None
            External forcing function
        n_frames : int
            Number of frames to output
        initial_perturbation : array or None
            Initial displacement from equilibrium

        Returns:
        --------
        t_frames : array
            Time points
        y_frames : array
            Shape (n_frames, n_nodes) - vertical positions at each frame
        """
        # Initial conditions
        if initial_perturbation is None:
            y0 = self.y_eq.copy()
        else:
            y0 = self.y_eq + initial_perturbation

        y_dot0 = np.zeros(self.n_nodes)

        # Boundary conditions
        y0[0] = self.y_eq[0]
        y0[-1] = self.y_eq[-1]

        state0 = np.concatenate([y0, y_dot0])

        # Time points for output
        t_eval = np.linspace(t_span[0], t_span[1], n_frames)

        # Solve ODE system
        print(f"Simulating dynamics from t={t_span[0]:.2f}s to t={t_span[1]:.2f}s...")
        sol = solve_ivp(
            fun=lambda t, y: self.equations_of_motion(t, y, slackliner_forces),
            t_span=t_span,
            y0=state0,
            method='RK45',
            t_eval=t_eval,
            rtol=1e-6,
            atol=1e-8,
        )

        print(f"Simulation complete. Generated {len(sol.t)} frames.")

        # Extract positions
        n = self.n_nodes
        y_frames = sol.y[:n, :].T  # Shape: (n_frames, n_nodes)

        return sol.t, y_frames


def create_slackliner_impulse(position, magnitude, duration=0.1):
    """
    Create an impulse force function for a slackliner jumping.

    Parameters:
    -----------
    position : float
        Horizontal position of the impulse
    magnitude : float
        Force magnitude (Newtons)
    duration : float
        Duration of the impulse (seconds)

    Returns:
    --------
    force_func : callable
        Function f(t, x, y) -> force array
    """
    def force_func(t, x, y):
        forces = np.zeros_like(x)
        if t < duration:
            # Find nearest node
            idx = np.argmin(np.abs(x - position))
            # Apply force with smooth ramp
            forces[idx] = magnitude * np.sin(np.pi * t / duration)
        return forces

    return force_func


def create_slackliner_oscillation(position, frequency, amplitude, phase=0):
    """
    Create an oscillating force function for a slackliner bouncing.

    Parameters:
    -----------
    position : float
        Horizontal position
    frequency : float
        Oscillation frequency (Hz)
    amplitude : float
        Force amplitude (Newtons)
    phase : float
        Phase offset (radians)

    Returns:
    --------
    force_func : callable
        Function f(t, x, y) -> force array
    """
    omega = 2 * np.pi * frequency

    def force_func(t, x, y):
        forces = np.zeros_like(x)
        # Find nearest node
        idx = np.argmin(np.abs(x - position))
        forces[idx] = amplitude * np.sin(omega * t + phase)
        return forces

    return force_func


def initial_pluck(position, displacement, width=2.0):
    """
    Create an initial pluck perturbation (Gaussian).

    Parameters:
    -----------
    position : float
        Horizontal position of pluck
    displacement : float
        Maximum vertical displacement (meters)
    width : float
        Width of the pluck (meters)

    Returns:
    --------
    perturbation_func : callable
        Function f(x) -> displacement array
    """
    def perturbation_func(x):
        return displacement * np.exp(-((x - position) / width)**2)

    return perturbation_func
