from .api import Constraints, Rig

def main():
    c = Constraints(gap_length=100, anchor_tension=1000)
    c.add_slackliner(position=50, mass=90)

    print({c.to_json()})
    c = Constraints.from_json(c.to_json())
    print({c.to_json()})

    r = c.rig()
    print(r.to_json())
    r = Rig.from_json(r.to_json())
    print(r.to_json())

