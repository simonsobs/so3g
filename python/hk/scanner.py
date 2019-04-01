import so3g
from spt3g import core
import numpy as np

class HKScanner:
    """Module that scans and reports on HK archive contents and compliance.

    """
    def __init__(self):
        self.session_id = None
        self.providers = {}
        self.stats = {
            'n_hk': 0,
            'n_other': 0,
            'n_session': 0,
        }

    def report(self):
        core.log_info('Report for session_id %i:\n' % self.session_id +
                      str(self.stats) + '\n' +
                      str(self.providers) + '\nEnd report.',
                      unit='HKScanner')
        self.session_id = None

    def __call__(self, f):
        """Processes a frame.  Only Housekeeping frames will be examined;
        other frames will simply be counted.  All frames are passed
        through unmodified.

        """
        if f.type == core.G3FrameType.EndProcessing:
            self.report()
            return [f]

        if f.type != core.G3FrameType.Housekeeping:
            self.stats['n_other'] += 1
            return f

        self.stats['n_hk'] += 1

        if f['hkagg_type'] == so3g.HKFrameType.session:
            session_id = f['session_id']
            if self.session_id is not None:
                if self.session_id != session_id:
                    self.report()
            if self.session_id is None:
                core.log_info('New HK Session id = %i, timestamp = %i' %
                              (session_id, f['start_time']), unit='HKScanner')
                self.session_id = session_id
                self.stats['n_session'] += 1

        elif f['hkagg_type'] == so3g.HKFrameType.status:
            # Have any providers disappeared?
            now_prov_id = [p['prov_id'].value for p in f['providers']]
            for p, info in self.providers.items():
                if p not in now_prov_id:
                    info['active'] = False
            
            # New providers?
            for p in now_prov_id:
                info = self.providers.get(p)
                if info is not None:
                    if not info['active']:
                        core.log_warn('prov_id %i came back to life.' % p,
                                      unit='HKScanner')
                        info['n_active'] += 1
                        info['active'] = True
                else:
                    self.providers[p] = {
                        'active': True, # Currently active (during processing).
                        'n_active': 1,  # Number of times this provider id became active.
                        'n_frames': 0,  # Number of data frames.
                        'timestamp_init': f['timestamp'],  # Timestamp of provider appearance
                        'timestamp_data': None, # Timestamp of most recent data frame.
                        'ticks': 0,   # Total number of timestamps in all blocks.
                        'span': None, # (earliest_time, latest_time)
                    }

        elif f['hkagg_type'] == so3g.HKFrameType.data:
            info = self.providers[f['prov_id']]
            info['n_frames'] += 1
            t_this = f['timestamp']
            if info['timestamp_data'] is None:
                t_ref = info['timestamp_init']
                if t_this < t_ref:
                    core.log_warn('data timestamp (%.1f) precedes provider '
                                  'timestamp by %f seconds.' % (t_this, t_this - t_ref),
                                  unit='HKScanner')
            elif t_this <= info['timestamp_data']:
                core.log_warn('data frame timestamps are not strictly ordered.',
                              unit='HKScanner')
            info['timestamp_data'] = t_this # update

            t_check = []
            for b in f['blocks']:
                if len(b.t):
                    if info['span'] is None:
                        info['span'] = b.t[0], b.t[-1]
                    else:
                        t0, t1 = info['span']
                        info['span'] = min(b.t[0], t0), max(b.t[-1], t1)
                    t_check.append(b.t[0])
                info['ticks'] += len(b.t)
            if len(t_check) and abs(min(t_check) - t_this) > 60:
                core.log_warn('data frame timestamp (%.1f) does not correspond to '
                              'data timestamp vectors (%.1f) .' % (t_this, t_check),
                              unit='HKScanner')
                
        else:
            core.log_warn('Weird hkagg_type: %i' % f['hkagg_type'],
                          unit='HKScanner')

        return [f]

if __name__ == '__main__':
    # Run me on a G3File containing a Housekeeping stream.
    core.set_log_level(core.G3LogLevel.LOG_INFO)
    import sys
    for f in sys.argv[1:]:
        p = core.G3Pipeline()
        p.Add(core.G3Reader(f))
        p.Add(HKScanner())
        p.Run()
