from collections import namedtuple

# A slackline parameter set is given by the following constants:
#   m: the mass per meter of the webbing
#   K: newtons per 100% stretch in webbing
#   g: the acceleration due to gravity
Slackline = namedtuple("Slackline", ["name", "m", "g", "K"])

# Leighton provided the following specs:
# TODO: proper references for these
dyneemite_pro = Slackline(name="Dyneemite Pro", m=0.088, g=9.81, K=2500*100.0)

def list_slacklines():
    return [dyneemite_pro]
