import so3g
from spt3g import core
import numpy as np

from so3g import hk

class HKScanner:
    """Module that scans and reports on HK archive contents and compliance.
    
    Attributes:
      stats (dict): A nested dictionary of statistics that are updated as
        frames are processed by the module.  Elements:

        - ``n_hk`` (int): The number of HK frames encountered.
        - ``n_other`` (int): The number of non-HK frames encountered.
        - ``n_session`` (int): The number of distinct HK sessions
          processed.
        - ``concerns`` (dict): The number of warning (key ``n_warning``)
          and error (key ``n_error``) events encountered.  The detail
          for such events is logged to ``spt3g.core.log_warning`` /
          ``log_error``.
    """
    def __init__(self):
        self.session_id = None
        self.providers = {}
        self.stats = {
            'n_hk': 0,
            'n_other': 0,
            'n_session': 0,
            'concerns': {
                'n_error': 0,
                'n_warning': 0
            },
            'versions': {}, # Number of frames counted for each hkagg_version.
        }

    def report_and_reset(self):
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
            self.report_and_reset()
            return [f]

        if f.type != core.G3FrameType.Housekeeping:
            self.stats['n_other'] += 1
            return f

        self.stats['n_hk'] += 1
        vers = f.get('hkagg_version', 0)
        self.stats['versions'][vers] = self.stats['versions'].get(vers, 0) + 1

        if f['hkagg_type'] == so3g.HKFrameType.session:
            session_id = f['session_id']
            if self.session_id is not None:
                if self.session_id != session_id:
                    self.report_and_reset()  # note this does clear self.session_id.
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
                        self.stats['concerns']['n_warning'] += 1
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
            vers = f.get('hkagg_version', 0)

            info['n_frames'] += 1
            t_this = f['timestamp']
            if info['timestamp_data'] is None:
                t_ref = info['timestamp_init']
                if t_this < t_ref:
                    core.log_warn('data timestamp (%.1f) precedes provider '
                                  'timestamp by %f seconds.' % (t_this, t_this - t_ref),
                                  unit='HKScanner')
                    self.stats['concerns']['n_warning'] += 1
            elif t_this <= info['timestamp_data']:
                core.log_warn('data frame timestamps are not strictly ordered.',
                              unit='HKScanner')
                self.stats['concerns']['n_warning'] += 1
            info['timestamp_data'] = t_this # update

            t_check = []

            if vers == 0:
                blocks = f['blocks']
                block_timef = lambda block: block.t
                block_itemf = lambda block: [(k, block.data[k]) for k in block.data.keys()]
            elif vers == 1:
                blocks = [f[k] for k in f.keys() if isinstance(f[k], so3g.IrregBlock)]
                block_timef = lambda block: np.array([t.time / core.G3Units.seconds for t in b.times])
                block_itemf = lambda block: [(k, block[k]) for k in block.keys()]

            for b in blocks:
                times = block_timef(b)
                if len(times):
                    if info['span'] is None:
                        info['span'] = times[0], times[-1]
                    else:
                        t0, t1 = info['span']
                        info['span'] = min(times[0], t0), max(times[-1], t1)
                    t_check.append(times[0])
                info['ticks'] += len(times)
                for k,v in block_itemf(b):
                    if len(v) != len(times):
                        core.log_error('Field "%s" has %i samples but .t has %i samples.' %

                                       (k, len(v), len(times)))
                        self.stats['concerns']['n_error'] += 1
            if len(t_check) and abs(min(t_check) - t_this) > 60:
                core.log_warn('data frame timestamp (%.1f) does not correspond to '
                              'data timestamp vectors (%s) .' % (t_this, t_check),
                              unit='HKScanner')
                self.stats['concerns']['n_warning'] += 1
                
        else:
            core.log_warn('Weird hkagg_type: %i' % f['hkagg_type'],
                          unit='HKScanner')
            self.stats['concerns']['n_warning'] += 1

        return [f]

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--translate', action='store_true')
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()

    # The report is displayed at level LOG_INFO.
    core.set_log_level(core.G3LogLevel.LOG_INFO)

    # Run me on a G3File containing a Housekeeping stream.
    for f in args.files:
        p = core.G3Pipeline()
        p.Add(core.G3Reader(f))
        if args.translate:
            p.Add(hk.HKTranslator())
        p.Add(HKScanner())
        p.Run()
