import so3g as ss
from spt3g import core

ss.TestClass().runme()
ss.greet()

print('Writing an HKInfo... maybe.')

w = core.G3Writer('out.g3')
f = core.G3Frame()
hk = ss.HKInfo()
hk.hk_source = 'hwp'
hk.session_id = 1324098
f['hk'] = hk
w.Process(f)

