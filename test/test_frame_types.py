import so3g as ss
from spt3g import core

ss.TestClass().runme()
ss.greet()

w = core.G3Writer('out.g3')
f = core.G3Frame()
f.type = core.G3FrameType.Scan
f['testv'] = core.G3VectorDouble([1., 2., 3., ])
w.Process(f)


