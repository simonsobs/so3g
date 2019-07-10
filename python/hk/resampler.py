import so3g
from spt3g import core
import numpy as np

class HKResampler:
    """Resample an HK data frame.

    Attributes
    ----------
    scheme: an object that has a method named `resample` that takes as an
        argument a frame block and returns the resampled block.
    """

    def __init__(self, scheme):
        self.scheme = scheme

    def __call__(self, f):
        """Processes a frame. Only Housekeeping frames with key `"hkagg_type"`
        set to `so3g.HKFrameType.data` will be resampled; any other frame will
        be returned as-is.

        Args:
            f: The frame to be processed

        Returns:
            The downsampled frame.
        """
        if f.type != core.G3FrameType.Housekeeping:
            return f
        if f["hkagg_type"] != so3g.HKFrameType.data:
            return f

        sess = so3g.hk.HKSessionHelper(f["session_id"])
        timestamp = min([b.t[0] for b in f["blocks"]])
        f_out = sess.data_frame(f["prov_id"], timestamp)
        for b in f["blocks"]:
            f_out["blocks"].append(self.scheme.resample(b))
        return [f_out]


class ResampleStep:
    """This resampling scheme simply reslices the input arrays with a step 
    size set by the user.

    Attributes:
        step: Positive, integer step size for resampling
    """
    def __init__(self, step):
        self.step = step

    def resample(self, block):
        """Resample an HK block.

        Parameters:
            block: A frame block.

        Returns:
            The resampled block.
        """
        ret = so3g.IrregBlockDouble()
        ret.t = np.array(block.t)[::self.step]
        for k in block.data.keys():
            ret.data[k] = np.array(block.data[k])[::self.step]
        return ret


class ResampleMedianMinMax:
    """This resampling scheme provides the median, minimum and maximum values in
    the specified window size.
    
    The median value is returned in a field with the same name; the others are
    the field name with `_min` and `_max` appended.

    For the timestamp, the median is returned.

    THIS CLASS IS STILL A DEMO: no attempt has been made to make it efficient
    and the min/max part hasn't been verified.

    Attributes:
        size: The length of each sample, in seconds.
    """
    def __init__(self, size):
        self.size = size

    def resample(self, block):
        ret = so3g.IrregBlockDouble()
        n = int((block.t[-1] - block.t[0]) / self.size) + 1
        ret.t = np.empty(n)
        for k in block.data.keys():
            ret.data[k] = np.empty(n)
            ret.data["%s_min" % k] = np.empty(n)
            ret.data["%s_max" % k] = np.empty(n)
        i_start = 0
        t_now = block.t[0] + self.size
        for j in range(n):
            i_end = i_start + 1
            if i_end < len(block.t):
                while block.t[i_end] < t_now:
                    if i_end + 1 < len(block.t):
                        i_end += 1
                    else:
                        break
            ret.t[j] = t_now - self.size / 2
            for k in block.data.keys():
                d = np.array(block.data[k])[i_start:i_end]
                if len(d):
                    ret.data[k][j] = np.median(d)
                    ret.data["%s_min" % k][j] = np.min(d)
                    ret.data["%s_max" % k][j] = np.max(d)
                else:
                    ret.data[k][j] = 0
                    ret.data["%s_min" % k][j] = 0
                    ret.data["%s_max" % k][j] = 0
            t_now += self.size
            i_start = i_end
        return ret


if __name__ == '__main__':
    import so3g
    import sys

    core.set_log_level(core.G3LogLevel.LOG_INFO)

    path = sys.argv[1:]
    p = core.G3Pipeline()
    p.Add(core.G3Reader(path))
    p.Add(so3g.hk.HKScanner())
    p.Add(so3g.hk.HKReframer(target=600.0))
    p.Add(so3g.hk.HKResampler(scheme=so3g.hk.ResampleMedianMinMax(3.0)))
    p.Add(so3g.hk.HKScanner())
    p.Add(core.G3Writer, filename="out_medianminmax.g3")
    p.Run()

    p = core.G3Pipeline()
    p.Add(core.G3Reader(path))
    p.Add(so3g.hk.HKScanner())
    p.Add(so3g.hk.HKReframer(target=600.0))
    p.Add(so3g.hk.HKResampler(scheme=so3g.hk.ResampleStep(10)))
    p.Add(so3g.hk.HKScanner())
    p.Add(core.G3Writer, filename="out_step.g3")
    p.Run()
