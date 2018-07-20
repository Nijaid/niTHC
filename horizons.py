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
                    if self.horizons.has_key(it):
                        continue
                    c = (float(dl_[i]) for i in [2,3,4])
                    d = (float(dl_[i+1]) - float(l_[i])
                            for i in [14, 16, 18])
                    self.horizons[it] = Ellipsoid(c, d)
