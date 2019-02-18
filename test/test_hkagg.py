import so3g
from spt3g import core
import numpy as np

test_file = 'hk_out.g3'

# Write a stream of HK frames.
# (Inspect the output with 'spt3g-dump hk_out.g3 so3g'.)
print('Streaming to %s' % test_file)
w = core.G3Writer(test_file)

# Create something to help us track the aggregator session.
hksess = so3g.hkagg.HKSession(session_id=1234,
                              description="Test HK data.")

# Register a data provider.
prov_id = hksess.add_provider(
    description='Fake data for the real world.')

# Start the stream -- write the initial session and status frames.
f = hksess.session_frame()
w.Process(f)
f = hksess.status_frame()
w.Process(f)

# Now make a data frame.
f = hksess.data_frame(prov_id=prov_id) 

# Add a data block.
hk = so3g.IrregBlockDouble()
hk.prefix = 'hwp_'
hk.data['position'] = [1,2,3,4,5]
hk.data['speed'] = [1.2, 1.2, 1.2, 1.2, 1.2]
hk.t = [0,1,2,3,4]
f['blocks'].append(hk)

w.Process(f)

del w
print('Stream closed.\n\n')

# Now play them back...
print('Reading back:')
for f in core.G3File(test_file):
    ht = f.get('hkagg_type')
    if ht == so3g.HKFrameType.session:
        print('Session: %i' % f['session_id'])
    elif ht == so3g.HKFrameType.status:
        print('  Status update: %i providers' % (len(f['providers'])))
    elif ht == so3g.HKFrameType.data:
        print('  Data: %i blocks' % len(f['blocks']))
        for block in f['blocks']:
            for k,v in block.data.items():
                print('    %s%s' % (block.prefix, k), v)
    
