#import sys
#sys.path.append('.')
import so3g as ss
from spt3g import core

r = ss.G3IndexedReader('/home/mhasse/data/spt3g/whitenoisesample.g3')
for i in range(10):
    pos = r.Tell()
    f = r.Process(None)[0]
    print(f.type)
    if f.type == core.G3FrameType.Wiring:
        print('Save')
        w_pos = pos
        
print ('Replay from wiring?')
r.Seek(w_pos)

for i in range(4):
    print(r.Process(None)[0].type)

