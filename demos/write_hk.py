# Note this demo is included directly from the package docs!
#
# This code generates a G3 file containing some telescope pointing
# data in the "SO HK" format.  When expanding it, check the SO HK
# format description to make sure your frame stream is compliant.

import time
import numpy as np

from spt3g import core
from so3g import hk

# Start a "Session" to help generate template frames.
session = hk.HKSessionHelper(hkagg_version=2)

# Create an output file and write the initial "session" frame.  If
# you break the data into multiple files, you write the session frame
# at the start of each file.
writer = core.G3Writer('hk_example.g3')
writer.Process(session.session_frame())

# Create a new data "provider".  This represents a single data
# source, sending data for some fixed list of fields.
prov_id = session.add_provider('platform')

# Whenever there is a change in the active "providers", write a
# "status" frame.
writer.Process(session.status_frame())

# Write, like, 10 half-scans.
frame_time = time.time()
v_az = 1.5  # degrees/second
dt = 0.001  # seconds
halfscan = 10  # degrees

for i in range(10):
    # Number of samples
    n = int(halfscan / v_az / dt)
    # Vector of unix timestamps
    t = frame_time + dt * np.arange(n)
    # Vector of az and el
    az = v_az * dt * np.arange(n)
    if i % 2:
      az = -az
    el = az * 0 + 50.

    # Construct a "block", which is a named G3TimesampleMap.
    block = core.G3TimesampleMap()
    block.times = core.G3VectorTime([core.G3Time(_t * core.G3Units.s)
                                     for _t in t])
    block['az'] = core.G3VectorDouble(az)
    block['el'] = core.G3VectorDouble(el)

    # Create an output data frame template associated with this
    # provider.
    frame = session.data_frame(prov_id)

    # Add the block and block name to the frame, and write it.
    frame['block_names'].append('pointing')
    frame['blocks'].append(block)
    writer.Process(frame)

    # For next iteration.
    frame_time += n * dt
