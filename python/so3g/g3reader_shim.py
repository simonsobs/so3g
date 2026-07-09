import warnings

from spt3g import core


class G3IndexedReader(core.G3Reader):
    def __init__(self, *a, **kw):
        warnings.warn("so3g.G3IndexedReader is deprecated and will be removed "
                      "in a future version; use spt3g.G3Reader (.seek/.tell).",
                      DeprecationWarning, stacklevel=2)
        return super().__init__(*a, **kw)
    def Seek(self, *a, **kw):
        return super().seek(*a, **kw)
    def Tell(self, *a, **kw):
        return super().tell(*a, **kw)

