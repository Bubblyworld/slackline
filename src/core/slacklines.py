from collections import namedtuple
from box import Box

# A slackline parameter set is given by the following constants:
#   m: the mass per meter of the webbing
#   K: newtons per 100% stretch in webbing
#   g: the acceleration due to gravity
class Slackline:
    def __init__(self, name, m, g, K):
        self.name = name
        self.m = m
        self.g = g
        self.K = K

    def __repr__(self):
        return "Slackline(name=%s, m=%d, g=%d, K=%d)" % (self.name, self.m, self.g, self.K)

    def __str__(self):
        return self.name

    def to_box(self):
        return Box({
            "name": self.name,
            "m": self.m,
            "g": self.g,
            "K": self.K,
        })

    @classmethod
    def from_box(cls, box):
        return cls(
            name=box.name,
            m=box.m,
            g=box.g,
            K=box.K,
        )

# Leighton provided the following specs:
# TODO: proper references for these
dyneemite_pro = Slackline(name="Dyneemite Pro", m=0.088, g=9.81, K=2500*100.0)

def list_slacklines():
    return [dyneemite_pro]
