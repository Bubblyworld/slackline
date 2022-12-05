from collections import namedtuple

# A slackline parameter set is given by the following constants:
#   m: the mass per meter of the webbing
#   K: newtons per 100% stretch in webbing
#   g: the acceleration due to gravity
Slackline = namedtuple("Slackline", ["m", "g", "K"])

# Leighton provided the following specs:
dyneemite_pro = Slackline(m=0.044, g=9.81, K=2500*100.0)
