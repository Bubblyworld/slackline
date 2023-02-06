"""
This module provides a high-level API for modelling slacklines. To use it,
first construct a Constaints object, which describes the physical properties
of the slackline and any slackliners that are on it. Then call the "rig" method
to obtain a collection of numpy data describing the length/tension/drop/angle
of the slackline at each point along the gap.

Constraints objects have modifiers that can be used to move slackliners, and
change the length of the line. These modifiers are designed to minimise the
amount of precomputation that happens in the background - use them to avoid
unnecessary work.
"""
from src.core.integrator import integrate, integrate_length_tension, integrate_natural_length
from src.core.calculus import first_order_euler_lagrange, mass_boundary_conditions
from src.core.slacklines import Slackline, list_slacklines, dyneemite_pro
from src.core.lagrangians import ideal
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
from box import Box
import json


class Rig:
    """
    A Rig object is a collection of numpy arrays that describe the length,
    tension, drop, and angle of the slackline at each point along the gap.
    It has the following attributes:
      - x: the horizontal distance along the gap in meters
      - n(x): the natural length of the slackline in meters
      - l(x): the arclength of the slackline in meters
      - y(x): the vertical drop of the slackline in meters
      - T(x): the tension in the slackline in Newtons
      - A(x): the angle of the slackline in degrees
    """
    def __init__(self, x, n, l, y, T, A):
        self.x = x
        self.n = n
        self.l = l
        self.y = y
        self.T = T
        self.A = A

    def to_dict(self):
        """
        Utility method that returns a dictionary representation of the rig.
        """
        return {
            "x": self.x.tolist(),
            "n": self.n.tolist(),
            "l": self.l.tolist(),
            "y": self.y.tolist(),
            "T": self.T.tolist(),
            "A": self.A.tolist(),
        }

    def to_json(self):
        """
        Utility method that returns a JSON representation of the rig.
        """
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json):
        """
        Utility method that constructs a Rig object from a JSON representation.
        """
        box = Box.from_json(json)
        return cls(
            x=np.array(box.x),
            n=np.array(box.n),
            l=np.array(box.l),
            y=np.array(box.y),
            T=np.array(box.T),
            A=np.array(box.A),
        )

    def plot(self):
        """
        Utility method that plots all of the rig's curves with matplotlib. You
        will need to call plt.show() to display the plot.
        """
        plt.figure(figsize=(12, 12))

        plt.subplot(221)
        plt.plot(self.x, self.y)
        plt.title("Curve")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")

        plt.subplot(222)
        plt.plot(self.x, self.n, label="n")
        plt.plot(self.x, self.l, label="l")
        plt.title("Natural and Arc Lengths")
        plt.xlabel("x (m)")
        plt.ylabel("length (m)")
        plt.legend()

        plt.subplot(223)
        plt.plot(self.x, self.T)
        plt.title("Tension")
        plt.xlabel("x (m)")
        plt.ylabel("T (N)")

        plt.subplot(224)
        plt.plot(self.x, self.A)
        plt.title("Angle")
        plt.xlabel("x (m)")
        plt.ylabel("A (deg)")

class Constraints:
    """
    A Constraints object describes the physical properties of a slackline and
    any slackliners that are on it.
    """
    def __init__(self, slackline=dyneemite_pro, gap_length=100, anchor_tension=1000):
        """
        Create a new Constraints object. Note that anchor_tension describes
        the *standing* anchor tension, i.e. before any slackliners are on the
        line.

        slackline: a Slackline object (default Dyneemite Pro)
        gap_length: the length of the gap in meters (default 100m)
        anchor_tension: standing anchor tension in newtons (default 1000N)
        """
        self.slackline = slackline
        self.gap_length = gap_length
        self.anchor_tension = anchor_tension
        self.slackliners = [] # list of (x, M) tuples

        # Cached precomputation:
        self.lagrangian = ideal
        self.foel = first_order_euler_lagrange(self.lagrangian)
        self.mbcs = [] # boundary value conditions for slackliners

    def add_slackliner(self, position, mass):
        """
        Adds a slackliner to the line. This will modify the Constraints object
        in-place, so that subsequent calls to rig() will include the slackliner.

        position: the horizontal position of the slackliner in meters
        mass: the mass of the slackliner in kilograms
        """
        if position < 0 or position > self.gap_length:
            raise ValueError("position must be between 0 and gap_length")
        if mass <= 0:
            raise ValueError("mass must be positive")

        self.slackliners.append((position, mass))
        self.mbcs.append(mass_boundary_conditions(
            self.lagrangian,
            position,
            self.slackline.m,
            self.slackline.g,
            self.slackline.K,
            mass,
        ))

    def rig(self):
        """
        Returns a Rig object that describes the rigged slackline at each point
        along the gap.
        """
        # First we compute the length-tension curve, as we know the gap length
        # and standing anchor tension. For this, we don't include any of the
        # slackliners, as you don't rig slacklines while people are standing on
        # them:
        x, y, n, y_x, n_x = integrate_length_tension(
            self.lagrangian,
            self.slackline,
            self.gap_length,
            self.anchor_tension,
            foel=self.foel,
        )

        if len(self.slackliners) > 0:
            # This gives us the natural length of the rigged slackline, which we
            # can use to compute the curves once the slackliners are on the line:
            natural_length = n[-1]
            x, y, n, y_x, n_x = integrate_natural_length(
                self.lagrangian,
                self.slackline,
                self.gap_length,
                natural_length,
                masses=self.slackliners,
                foel=self.foel,
                mbcs=self.mbcs,
            )

        # The arclength can be computed differentially using Pythagoras' law:
        dx = x[1] - x[0] # assuming constant step size!
        dl = np.sqrt(1 + y_x**2) * dx
        l = np.cumsum(dl)

        # The tension in the line can be computed as follows:
        dn = n_x * dx
        T = self.slackline.K * (dl/dn - 1)

        # Finally, the angle is related to the vertical drop:
        A = np.abs(np.arctan(y_x)) * 180 / np.pi

        return Rig(x, n, l, y, T, A)

    def to_json(self):
        """
        Utility method that returns a JSON representation of the Constraints
        object.
        """
        return Box({
            "slackline": self.slackline.to_box(),
            "gap_length": self.gap_length,
            "anchor_tension": self.anchor_tension,
            "slackliners": self.slackliners,
        }).to_json()

    @classmethod
    def from_json(cls, json):
        """
        Utility method that constructs a Constraints object from a JSON
        representation.
        """
        if type(json) == str:
            json = json.loads(json)
        box = Box(json)
        constraints = cls(
            slackline=Slackline.from_box(box.slackline),
            gap_length=box.gap_length,
            anchor_tension=box.anchor_tension,
        )
        for position, mass in box.slackliners:
            constraints.add_slackliner(position, mass)
        return constraints
