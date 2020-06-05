import so3g
from spt3g import core
import numpy as np
import pickle
import datetime, time
import sys, os
import warnings
import argparse

def g3_to_array(g3file, verbose=False):
    """
    Takes a G3 file output from the SMuRF archiver and reads to a numpy array.

    Parameters
    ----------
    g3file : full path to the G3 file
    verbose : OPTIONAL choice for verbose output (-v, --verbosity)

    Returns
    -------
    times : array of G3Time objects
    data : array of arrays, where each internal array is a SMuRF channel
    """
    frames = [fr for fr in core.G3File(g3file)]
    
    data=[]

    frametimes = []
    for frame in frames:
        if frame.type == core.G3FrameType.Scan:
            frametime = frame['data'].times()
            frametimes.append(frametime)

    if frametimes == []:
        warnings.warn('Wrong frame type')

    strtimes = np.hstack(frametimes)
    
    times = []
    for strt in strtimes:
        t=core.G3Time(strt).time/core.G3Units.s
        times.append(t)
    times = np.asarray(times)
    
    channums = []

    i=0
    while i<len(frames):
        if verbose:
            print('Trying frame %i'%i)
        frametype = frames[i].type
        if frametype == core.G3FrameType.Scan:
            for chan in frames[i]['data'].keys():
                channums.append(int(chan.strip('r')))
            break
        else:
            i+=1
    if verbose:
        print('Channel numbers obtained')

    channums.sort()
    for ch in channums:
        if verbose:
            print('Adding channel %s'%ch)
        chdata = []
        for frame in frames:
            if frame.type == core.G3FrameType.Scan:
                framedata = frame['data']['r'+format(ch, "04")]
                chdata.append(framedata)
        chdata_all = np.hstack(chdata)
        data.append(chdata_all)

    biases = []
    biasnums = []
    for num in frames[i]['tes_biases'].keys():
        biasnums.append(num)
    biasnums.sort()
    for b in biasnums:
        if verbose:
            print('Adding bias {}'.format(b))
        bias = []
        for frame in frames:
            if frame.type == core.G3FrameType.Scan:
                biasdata = frame['tes_biases'][b]
                bias.append(biasdata)
        bias_all = np.hstack(bias)
        biases.append(bias_all)
    biases = np.asarray(biases)
    data = np.asarray(data)
    return times, data, biases


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('g3file', type=str, metavar='G3', help='full path to g3 file')
    parser.add_argument('-v', '--verbosity', action='store_true', default=False)

    args = parser.parse_args()

    times, data, biases = g3_to_array(args.g3file, args.verbosity)
    newfile_name = args.g3file.split('/')[-1].strip('.g3')
    with open(newfile_name+'.pkl','wb') as f:
        pickle.dump({'times':times, 'data':data, 'tes_biases':biases},f)
    f.close()

