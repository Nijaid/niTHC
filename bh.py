from scivis.utils import Ellipsoid
from scidata.utils import ilocate

class BHHorizons(object):
    def __init__(self, root="./"):
        self.horizons = {}
        for gp in ilocate("BH_diagnostics.ah1.gp", root=root,
                            followlinks=True):
            with open(gp, "r") as f_:
                for line in f_:
                    if line[0]=='#':
                        continue
                    l_ = line.split()
                    it = int(l_[0])
                    if it in self.horizons:
                        continue
                    c = (float(l_[i]) for i in [2,3,4])
                    d = (float(l_[i+1]) - float(l_[i])
                            for i in [14, 16, 18])
                    self.horizons[it] = Ellipsoid(c, d)

    def horizon(self, plane, it):
        axis = {"xy": 2, "xz": 1, "yz": 0}[plane]
        if it in self.horizons.keys():
            h = self.horizons[it]
            return(h.slice(s=0., axis=axis))
        else:
            return(None)
