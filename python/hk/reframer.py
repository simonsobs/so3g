import so3g
from spt3g import core
import numpy as np

class _HKBlockBundle:
    def __init__(self):
        self.t = None

    def add(self, b):
        """Cache the data from IrregBlockDouble b."""
        if self.t is None:
            self.t = []
            self.chans = {k: [] for k in b.data.keys()}
        self.t.extend(list(b.t))
        for k in b.data.keys():
            self.chans[k].extend(b.data[k])

    def ready(self, flush_time):
        """Returns true if the current block has crossed flush_time."""
        return len(self.t) > 0 and self.t[-1] >= flush_time

    def rebundle(self, flush_time):
        """Return the IrregBlockDouble with all samples timestamped up to but
        not including flush_time.  Pass None to force rebundle all.

        If there are no samples to bundle, None is returned.
        """
        if len(self.t) == 0:
            return None
        if flush_time is None:
            idx = len(self.t)
        else:
            idx = 0
            while idx < len(self.t) and self.t[idx] < flush_time:
                idx += 1
        out = so3g.IrregBlockDouble()
        out.t = np.array(self.t[:idx])
        self.t = self.t[idx:]
        for k in self.chans.keys():
            out.data[k] = np.array(self.chans[k][:idx])
            self.chans[k] = self.chans[k][idx:]
        return out


class _HKProvBundle:
    def __init__(self, t0, dt):
        self.blocks = []
        self.field_map = {}
        self.sess = None
        self.prov_id = None
        self.t0 = t0
        self.dt = dt

    def add(self, f):
        if self.sess is None:
            self.sess = so3g.hk.HKSessionHelper(f['session_id'])
            self.prov_id = f['prov_id']
        for b in f['blocks']:
            chans = b.data.keys()
            if len(chans) == 0 or len(b.t) == 0:
                continue
            idx = self.field_map.get(chans[0], None)
            if idx is None:
                idx = len(self.blocks)
                self.blocks.append(_HKBlockBundle())
                for k in chans:
                    self.field_map[k] = idx
            self.blocks[idx].add(b)

    def ready(self):
        return all(b.ready(self.t0 + self.dt) for b in self.blocks)

    def rebundle(self, force=False):
        output = []
        while force or self.ready():
            if force:
                flush_time = None
            else:
                flush_time = self.t0 + self.dt
                self.t0 += self.dt
            blocks_out = [b.rebundle(flush_time) for b in self.blocks]
            blocks_out = [b for b in blocks_out if b is not None and len(b.t)>0]
            if len(blocks_out) == 0:
                break
            timestamp = min([b.t[0] for b in blocks_out])
            f = self.sess.data_frame(self.prov_id, timestamp)
            for b in blocks_out:
                f['blocks'].append(b)
            output.append(f)
        return output


class HKReframer:
    """A module to rebundle SO HK frames, to decrease or increase the
    number of samples in each data frame.

    """

    def __init__(self, target=60.):
        """
        Arguments:
          target (float): The target frame duration, in seconds.
        """
        self.target = target
        self.session_id = None
        self.providers = {}

    def flush(self, prov_id=None):
        """Flush all buffers; empty the provider list; reset session_id.
        Returns a list of flushed output frames.

        If prov_id is specified, then only that provider is flushed
        and popped from provider list.

        """
        if prov_id is not None:
            return self.providers.pop(prov_id).rebundle(True)
        core.log_info('Flushing session id %i' % (self.session_id),
                      unit='HKReframer')
        output = []
        for p in self.providers.values():
            output += p.rebundle(True)
        self.providers = {}
        self.session_id = None
        return output

    def __call__(self, f):
        """Processes a frame.  Only Housekeeping frames will be manipulated;
        others will be passed through untouched.

        """
        if f.type == core.G3FrameType.EndProcessing:
            return self.flush() + [f]

        if f.type != core.G3FrameType.Housekeeping:
            return f

        output = []

        if f['hkagg_type'] == so3g.HKFrameType.session:
            session_id = f['session_id']
            if self.session_id is not None:
                if self.session_id != session_id:
                    output += self.flush()  # this clears self.session._id.
                    output.append(f)
                else:
                    pass # Don't re-emit an on-going session frame.
            if self.session_id is None:
                core.log_info('New HK Session id = %i, timestamp = %i' %
                              (session_id, f['start_time']), unit='HKReframer')
                self.session_id = session_id
                output.append(f)

        elif f['hkagg_type'] == so3g.HKFrameType.status:
            # Only issue status if something has changed.
            changes = False
            # Flush any providers that are now expired.
            now_prov_id = [p['prov_id'].value for p in f['providers']]
            for p in list(self.providers.keys()):
                if p not in now_prov_id:
                    output += self.flush(p)
                    changes = True
            # Create bundlers for any new providers.
            for p in now_prov_id:
                if p not in self.providers:
                    t0 = f['timestamp']
                    self.providers[p] = _HKProvBundle(t0, self.target)
                    changes = True
            if changes:
                output.append(f)

        elif f['hkagg_type'] == so3g.HKFrameType.data:
            fb = self.providers[f['prov_id']]
            fb.add(f)
            if fb.ready():
                output += fb.rebundle()

        else:
            raise ValueError('Invalid hkagg_type')

        return output


if __name__ == '__main__':
    from so3g.hk import HKScanner, HKReframer
    import sys

    core.set_log_level(core.G3LogLevel.LOG_INFO)

    files = sys.argv[1:]
    p = core.G3Pipeline()
    p.Add(core.G3Reader(files))
    p.Add(HKScanner())
    rf = HKReframer()
    p.Add(rf)
    p.Add(HKScanner())
    p.Add(core.G3Writer, filename='out.g3')
    p.Run()
