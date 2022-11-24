import matplotlib.pyplot as plt
import numpy as np

# Constants:
g = 9.81 # m / s^2 (gravitational acceleration)
K = 2500 # N (newtons per 1% stretch in rope)
m = 0.044 # kg / m (mass density of rope in natural units)

def _curve(T_x, T_y, N):
    """
    The problem of determining the webbing's curve is easiest to solve in
    natural units, where we parameterise the curve by its untensioned length.
    Since we allow tension to vary, and the webbing stretches depending on
    tension, converting back into SI units is a bit tricky, and requires us to
    solve for the curve/tension before we can get anywhere.

    The following equations describes the curve u in terms of natural units n:
        a = (T_y - m * g * n) / T_x
        T = T_x * sqrt(1 + a^2)
        dx/dn = (T/100/K + 1) * T_x/T
        du/dn = a * dx/dn

    This function solves for the curve u(n) and the SI units x(n), where n
    \in [0, N].
    """
    # The natural units 'n':
    n = np.linspace(0, N, 1000)
    dn = n[1] - n[0]

    # The auxilliary variable 'a' and the tension 'T':
    a = (T_y - 2 * m * g * n) / T_x
    T = T_x * np.sqrt(1 + a**2)

    # The SI unit 'x':
    dx = (T/100/K + 1) * T_x/T * dn
    x = np.cumsum(dx)

    # The actual curve 'u':
    du = a * dx
    u = np.cumsum(du)

    return x, u, {"dx":dx, "du":du, "a":a, "T":T, "T_x":T_x, "T_y":T_y, "dn":dn}

def curve(T, D):
    """
    Solves for the curve of the webbing given T (the magnitude of the anchor
    tension in Newtons) and D (the length of the gap in meters).
    """
    def _angle(rads):
        """
        Finds the curve for a given angle of the anchor webbing in radians. As
        we don't know the length of the curve in natural units, we binary
        search on possible Ns until u(N)=0 to 3 decimal places.
        """
        # Binary search for the natural units 'N':
        T_x = T * np.cos(rads)
        T_y = T * np.sin(rads)
        N_min = 1e-6
        N_max = 1e6
        while N_max - N_min > 1e-3:
            N = (N_min + N_max) / 2
            x, u, md = _curve(T_x, T_y, N)
            if u[-1] > 0:
                N_min = N
            else:
                N_max = N

        # Return the curve:
        return x, u, md

    # We need to binary search for the correct angle in radians, since we don't
    # know which angle will give us an SI length (x) that is equal to D:
    rads_min = np.radians(1e-6)
    rads_max = np.radians(90)
    while True:
        rads = (rads_min + rads_max) / 2
        x, u, md = _angle(rads)
        if x[-1] < D:
            rads_min = rads
        else:
            rads_max = rads
        if abs(x[-1] - D) < 1e-3:
            break
    return x, u, md

def total_force(x, u, md):
    """
    Calculates the magnitude of total force on each point of the webbing.
    """
    a = md["a"]
    T = md["T"]
    T_x = md["T_x"]
    T_y = md["T_y"]
    dn = md["dn"]
    dx = md["dx"]
    du = md["du"]
    dl = np.sqrt(dx**2 + du**2)
    sin_theta = du / dl
    cos_theta = dx / dl

    f_x = T[1:] * cos_theta[1:] - T[:-1] * cos_theta[:-1]
    f_y = T[1:] * sin_theta[1:] - T[:-1] * sin_theta[:-1] + 2 * m * g * dn
    return np.sqrt(f_x**2 + f_y**2)

def main():
    D = 100 # m (length of the gap)
    Ts = [500, 1000, 1500, 2000] # N (tension in the anchor webbing)

    # Plot the curves:
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.set_xlabel("x (m)")
    ax1.set_ylabel("u (m)")
    ax1.set_title("Webbing Curve")

    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("f (N)")
    ax2.set_ylim(-0.1, 1)
    ax2.set_title("Total Force")

    ax3.set_xlabel("x (m)")
    ax3.set_ylabel("du/dx")

    ax4.set_xlabel("x (m)")
    ax4.set_ylabel("T (N)")
    ax4.set_title("Tension")

    # Plot the curve and tension for each T:
    for T in Ts:
        x, u, md = curve(T, D)
        f = total_force(x, u, md)
        ax1.plot(x, -u, label=f"T={T} N")
        ax2.plot(x[:-1], f, label=f"T={T} N")
        ax3.plot(x, md["du"]/md["dx"], label=f"T={T} N")
        ax4.plot(x, md["T"], label=f"T={T} N")

    # Adjust the axes so they don't overlap:
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    fig.tight_layout()
    plt.show()
