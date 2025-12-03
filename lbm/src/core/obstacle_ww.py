### ************************************************
### Class defining an obstacle in the lattice
class obstacle:
    def __init__(self, name):

        self.name   = name

    def fill(self, obs, boundary, area):

        self.obs     = obs
        self.boundary = boundary
        self.area     = area
