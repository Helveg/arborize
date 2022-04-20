import glia as g
from patch import p

class Synapse:

    def __init__(self, cell, section, point_process_name, attributes = {}, variant=None, type=None, source=None):
        self._cell = cell
        self._type = type
        self._section = section
        self._point_process_name = point_process_name
        self.source = source
        with g.context(pkg=cell.__class__.glia_package):
            if "gap" in point_process_name.lower():
                self._point_process_glia_name = g.resolve(point_process_name, variant=variant)
                self._point_process = g.insert(section, point_process_name, variant=variant)
            else:
                self._point_process_glia_name = "ExpSyn"
                self._point_process = p.ExpSyn(section(0.5))
        section.__ref__(self)
        for key, value in attributes.items():
            setattr(self._point_process, key, value)


    def __neuron__(self):
        return self._point_process.__neuron__()

    def stimulate(self, *args, **kwargs):
        return self._point_process.stimulate(*args, **kwargs)

    def record(self):
        return p.record(self._point_process._ref_i)

    def presynaptic(self, section, x=0.5, **kwargs):
        if self.source is None:
            nc = p.NetCon(section(x)._ref_v, self._point_process, sec=section, **kwargs)
            w = -1 if "GABA" in self._point_process_name else 1
            nc.weight[0] = 0.00001 * w
            return nc
        else:
            setattr(self._point_process, f"_ref_{self.source}", section(x)._ref_v)
