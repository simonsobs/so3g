import so3g
from spt3g import core
import time

class HKSessionHelper:
    """
    Helper class to produce G3Frame templates for creating streams of
    generic HK data.
    """
    def __init__(self, session_id=0, start_time=None,
                 description='No description provided.'):
        if start_time is None:
            start_time = time.time()
        self.session_id = session_id
        self.start_time = start_time
        self.description = description
        self.provs = {}
        self.next_prov_id = 0

    def add_provider(self, description='No provider description... provided'):
        prov_id = self.next_prov_id
        self.next_prov_id += 1
        self.provs[prov_id] = {'description': description}
        return prov_id

    def remove_provider(self, prov_id):
        del self.provs[prov_id]

    """
    Frame generators.
    """
        
    def session_frame(self):
        """
        Return the Session frame.  No additional information needs to be
        added to this frame.  This frame initializes the HK stream (so
        it should be the first frame of each HK file; it should
        precede all other HK frames in a network source).
        """
        f = core.G3Frame()
        f.type = core.G3FrameType.Housekeeping
        f['hkagg_type'] = so3g.HKFrameType.session
        f['session_id'] = self.session_id
        f['start_time'] = self.start_time
        f['description'] = self.description
        return f

    def status_frame(self, timestamp=None):
        """
        Return a Status frame template.  Before processing, the session
        manager should update with information about what providers
        are currently connected.
        """
        if timestamp is None:
            timestamp = time.time()
        f = core.G3Frame()
        f.type = core.G3FrameType.Housekeeping
        f['hkagg_type'] = so3g.HKFrameType.status
        f['session_id'] = self.session_id
        f['timestamp'] = timestamp
        provs = core.G3VectorFrameObject()
        for prov_id in sorted(self.provs.keys()):
            prov = core.G3MapFrameObject()
            prov['prov_id'] = core.G3Int(prov_id)
            prov['description'] = core.G3String(
                self.provs[prov_id]['description'])
            provs.append(prov)
        f['providers'] = provs
        return f
    
    def data_frame(self, prov_id, timestamp=None):
        """
        Return a Data frame template.  The prov_id must match the prov_id
        in one of the Provider blocks in the preceding status frame.
        The session manager should create and add IrregBlockDouble items
        to the 'blocks' list.
        """
        if timestamp is None:
            timestamp = time.time()
        f = core.G3Frame()
        f.type = core.G3FrameType.Housekeeping
        f['hkagg_type'] = so3g.HKFrameType.data
        f['session_id'] = self.session_id
        f['prov_id'] = prov_id
        f['timestamp'] = timestamp
        f['blocks'] = core.G3VectorFrameObject()
        return f
