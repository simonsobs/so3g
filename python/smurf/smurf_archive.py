import sqlalchemy as db
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref

from spt3g import core
import so3g
import datetime as dt
import os
from tqdm import tqdm
import numpy as np
import yaml
import ast
from collections import namedtuple
from enum import Enum


Base = declarative_base()
Session = sessionmaker()
num_bias_lines = 16

class Files(Base):
    """Table to store file indexing info"""
    __tablename__ = 'files'
    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String, nullable=False)
    start = db.Column(db.DateTime)
    stop = db.Column(db.DateTime)
    frames = db.Column(db.Integer)
    channels = db.Column(db.Integer)

type_key = ['Observation', 'Wiring', 'Scan']
class FrameType(Base):
    """Enum table for storing frame types"""
    __tablename__ = "frame_type"
    id = db.Column(db.Integer, primary_key=True)
    type_name = db.Column(db.String, unique=True, nullable=True)

class Frame(Base):
    """Table to store frame indexing info"""
    __tablename__ = 'frames'
    __table_args__ = (
        db.UniqueConstraint('file_id', 'frame_idx', name='_frame_loc'),
    )

    id = db.Column(db.Integer, primary_key=True)

    file_id = db.Column(db.Integer, db.ForeignKey('files.id'))
    file = relationship("Files")

    frame_idx = db.Column(db.Integer, nullable=False)
    offset = db.Column(db.Integer, nullable=False)

    type_name = db.Column(db.String, db.ForeignKey('frame_type.type_name'))
    frame_type = relationship('FrameType')

    time = db.Column(db.DateTime, nullable=False)

    # Specific to data frames
    samples = db.Column(db.Integer)
    channels = db.Column(db.Integer)
    start = db.Column(db.DateTime)
    stop = db.Column(db.DateTime)

    def __repr__(self):
        return f"Frame({self.type_name})<{self.location}>"


class TimingParadigm(Enum):
    G3Timestream = 1
    SmurfUnixTime = 2
    TimingSystem = 3
    Mixed = 4

def get_sample_timestamps(frame):
    """
    Gets timestamps of samples in a G3Frame. This will try to get highest
    precision first and move to lower precision methods if that fails.

    Args
    ------
        frame (core.G3Frame):
            A G3Frame(Scan) containing streamed detector data.

    Returns
    ---------
        times (np.ndarray):
            numpy array containing timestamps in seconds
        paradigm (TimingParadigm):
            Paradigm used to calculate timestamps.
    """
    if 'primary' in frame.keys():
        if False:
            # Do high precision timing calculation here when we have real data
            pass
        else:
            # Try to calculate the timestamp based on the SmurfProcessor's
            # "UnixTime" and the G3Timestream start time.  "UnixTime" is a
            # 32-bit nanosecond clock that steadily increases mod 2**32.
            unix_times = np.array(frame['primary']['UnixTime'])
            for i in np.where(np.diff(unix_times) < 0)[0]:
                # This corrects for any wrap around
                unix_times[i+1:] += 2**32
            times = frame['data'].start.time / core.G3Units.s \
                + (unix_times - unix_times[0]) / 1e9

            return times, TimingParadigm.SmurfUnixTime
    else:
        # Calculate timestamp based on G3Timestream.times(). Note that this
        # only uses the timestream start and end time, and assumes samples are
        # equispaced.
        times = np.array([t.time / core.G3Units.s
                          for t in frame['data'].times()])
        return times, TimingParadigm.G3Timestream

class SmurfArchive:
    def __init__(self, archive_path, db_path=None, echo=False):
        """
        Class to manage a smurf data archive.

        Args
        -----
            archive_path (path):
                Path to the data directory
            db_path (path, optional):
                Path to the sqlite file. Defaults to ``<archive_path>/frames.db``
            echo (bool, optional):
                If true, all sql statements will print to stdout.
        """
        if db_path is None:
            db_path = os.path.join(archive_path, 'frames.db')
        self.archive_path = archive_path
        self.engine = db.create_engine(f"sqlite:///{db_path}", echo=echo)
        Session.configure(bind=self.engine)
        self.Session = sessionmaker(bind=self.engine)
        Base.metadata.create_all(self.engine)

        # Defines frame_types
        self._create_frame_types()

    def _create_frame_types(self):
        session = self.Session()
        if not session.query(FrameType).all():
            print("Creating FrameType table...")
            for k in type_key:
                ft = FrameType(type_name=k)
                session.add(ft)
            session.commit()

    def add_file(self, path, session):
        """
        Indexes a single file and adds it to the sqlite database.

        Args
        ----
            path (path): Path of the file to index
        """

        frame_types = {
            ft.type_name: ft for ft in session.query(FrameType).all()
        }

        db_file = Files(path=path)
        session.add(db_file)

        reader = so3g.G3Reader(path)

        total_channels = 0
        file_start, file_stop = None, None
        frame_idx = 0
        while True:
            db_frame = Frame(frame_idx=frame_idx, file=db_file)
            db_frame.offset = reader.tell()

            frames = reader.Process(None)
            if not frames:
                break
            frame = frames[0]
            frame_idx += 1

            if str(frame.type) not in type_key:
                continue

            db_frame.frame_type = frame_types[str(frame.type)]

            timestamp = frame['time'].time / core.G3Units.s
            db_frame.time = dt.datetime.fromtimestamp(timestamp)

            data = frame.get('data')
            if data is not None:
                db_frame.samples = data.n_samples
                db_frame.channels = len(data)
                db_frame.start = dt.datetime.fromtimestamp(data.start.time / core.G3Units.s)
                db_frame.stop = dt.datetime.fromtimestamp(data.stop.time / core.G3Units.s)

                if file_start is None:
                    file_start = db_frame.start
                file_stop = db_frame.stop
                total_channels = max(total_channels, db_frame.channels)

            session.add(db_frame)

        db_file.start = file_start
        db_file.stop = file_stop
        db_file.channels = total_channels
        db_file.frames = frame_idx


    def index_archive(self, verbose=False, stop_at_error=False):
        """
        Adds all files from an archive to the sqlite database.

        Args
        ----
        verbose: bool
            Verbose mode
        stop_at_error: bool
            If True, will stop if there is an error indexing a file.
        """
        session = self.Session()
        indexed_files = [f[0] for f in session.query(Files.path).all()]

        files = []
        for root, _, fs in os.walk(self.archive_path):
            for f in fs:
                path = os.path.join(root, f)
                if path.endswith('.g3') and path not in indexed_files:
                    files.append(path)

        if verbose:
            print(f"Indexing {len(files)} files...")

        for f in tqdm(files):
            try:
                self.add_file(os.path.join(root, f), session)
                session.commit()
            except IntegrityError as e:
                # Database Integrity Errors, such as duplicate entries
                session.rollback()
                print(e)
            except RuntimeError as e:
                # End of stream errors, for G3Files that were not fully flushed
                session.rollback()
                print(f"Failed on file {f} due to end of stream error!")
            except Exception as e:
                # This will catch generic errors such as attempting to load
                # out-of-date files that do not have the required frame
                # structure specified in the TOD2MAPS docs.
                session.rollback()
                if stop_at_error:
                    raise e
                elif verbose:
                    print(f"Failed on file {f}:\n{e}")

        session.close()

    def load_data(self, start, end, show_pb=True, load_biases=True):
        """
        Loads smurf G3 data for a given time range. For the specified time range
        this will return a chunk of data that includes that time range.

        Args
        -----
            start (timestamp): start timestamp
            end   (timestamp): end timestamp
            show_pb (bool, optional): If True, will show progress bar.
            load_biases (bool, optional): If True, will return biases.

        Returns
        --------
            Returns a tuple ``SmurfData(times, data, primary, status, biases, timing_paradigm)``
            with the following fields:

                times (np.ndarray[samples]):
                    Array of unix timestamps for loaded data
                data (np.ndarray[channels, samples]):
                    Array of the squid phase in units of radians for each
                    channel with data in the specified time range. The index of
                    the array is the readout channel number.
                primary (Dict[np.ndarray]):
                    Dict of numpy arrays holding the "primary" data included in
                    the packet headers, including 'AveragingResetBits',
                    'Counter0', 'Counter1', 'Counter2', 'FluxRampIncrement',
                    'FluxRampOffset', 'FrameCounter', 'TESRelaySetting',
                    'UnixTime'
                status (SmurfStatus):
                    SmurfStatus object containing metadata info at the time of
                    the first Scan frame in the requested interval. If there
                    are no Scan frames in the interval, this will be None.
                biases (optional, np.ndarray[NTES, samples]):
                    An array containing the TES bias values.
                    If ``load_biases`` is False, this will be None.
                timing_paradigm(TimingParadigm):
                    Tells you the method used to extract timestamps from the
                    frame data.
        """
        session = self.Session()

        frames = session.query(Frame).filter(
            Frame.type_name == 'Scan',
            Frame.stop >= dt.datetime.fromtimestamp(start),
            Frame.start < dt.datetime.fromtimestamp(end)
        ).order_by(Frame.time)
        session.close()

        samples, channels = 0, 0
        num_frames = 0
        for f in frames:
            num_frames += 1
            samples += f.samples
            channels = max(f.channels, channels)

        timestamps = np.full((samples,), np.nan, dtype=np.float64)
        data = np.full((channels, samples), 0, dtype=np.int32)
        if load_biases:
            biases = np.full((num_bias_lines, samples), 0, dtype=np.int32)
        else:
            biases = None

        primary = {}

        cur_sample = 0
        cur_file = None
        timing_paradigm = None
        for frame_info in tqdm(frames, total=num_frames, disable=(not show_pb)):
            file = frame_info.file.path
            if file != cur_file:
                reader = so3g.G3Reader(file)
                cur_file = file

            reader.seek(frame_info.offset)
            frame = reader.Process(None)[0]
            nsamp = frame['data'].n_samples

            key_order = [int(k[1:]) for k in frame['data'].keys()]
            data[key_order, cur_sample:cur_sample + nsamp] = frame['data']

            if load_biases:
                bias_order = [int(k[-2:]) for k in frame['tes_biases'].keys()]
                biases[bias_order, cur_sample:cur_sample + nsamp] = frame['tes_biases']

            # Loads primary data
            if 'primary' in frame.keys():
                for k, v in frame['primary'].items():
                    if k not in primary:
                        primary[k] = np.zeros(samples, dtype=np.int64)
                    primary[k][cur_sample:cur_sample + nsamp] = v

            ts, paradigm = get_sample_timestamps(frame)
            if timing_paradigm is None:
                timing_paradigm = paradigm
            elif timing_paradigm != paradigm:
                timing_paradigm = TimingParadigm.Mixed

            timestamps[cur_sample:cur_sample + nsamp] = ts

            cur_sample += nsamp

        # Conversion from DAC counts to squid phase
        rad_per_count = np.pi / 2**15
        data = data * rad_per_count

        if len(timestamps) > 0:
            status = self.load_status(timestamps[0])
        else:
            status = None

        SmurfData = namedtuple('SmurfData', 'times data primary status biases timing_paradigm')
        if load_biases:
            return SmurfData(timestamps, data, primary, status, biases, timing_paradigm)
        else:
            return SmurfData(timestamps, data, primary, status, None, timing_paradigm)

    def load_status(self, time, show_pb=False):
        """
        Returns the status dict at specified unix timestamp.
        Loads all status frames between session start frame and specified time.

        Args:
            time (timestamp): Time at which you want the rogue status

        Returns:
            status (dict): Dictionary of rogue variables at specified time.
        """
        session = self.Session()
        session_start,  = session.query(Frame.time).filter(
            Frame.type_name == 'Observation',
            Frame.time <= dt.datetime.fromtimestamp(time)
        ).order_by(Frame.time.desc()).first()

        status_frames = session.query(Frame).filter(
            Frame.type_name == 'Wiring',
            Frame.time >= session_start,
            Frame.time <= dt.datetime.fromtimestamp(time)
        ).order_by(Frame.time)

        status = {}
        cur_file = None
        for frame_info in tqdm(status_frames.all(), disable=(not show_pb)):
            file = frame_info.file.path
            if file != cur_file:
                reader = so3g.G3Reader(file)
                cur_file = file
            reader.seek(frame_info.offset)
            frame = reader.Process(None)[0]
            status.update(yaml.safe_load(frame['status']))

        return SmurfStatus(status)


class SmurfStatus:
    """
    This is a class that attempts to extract essential information from the
    SMuRF status dictionary so it is more easily accessible. If the necessary
    information for an attribute is not present in the dictionary, the
    attribute will be set to None.

    Args
    -----
        status  : dict
            A SMuRF status dictionary

    Attributes
    ------------
        status : dict
            Full smurf status dictionary
        num_chans: int
            Number of channels that are streaming
        mask : Optional[np.ndarray]
            Array with length ``num_chans`` that describes the mapping
            of readout channel to absolute smurf channel.
        mask_inv : np.ndarray
            Array with dimensions (NUM_BANDS, CHANS_PER_BAND) where
            ``mask_inv[band, chan]`` tells you the readout channel for a given
            band, channel combination.
        freq_map : Optional[np.ndarray]
            An array of size (NUM_BANDS, CHANS_PER_BAND) that has the mapping
            from (band, channel) to resonator frequency. If the mapping is not
            present in the status dict, the array will full of np.nan.
        filter_a : Optional[np.ndarray]
            The A parameter of the readout filter.
        filter_b : Optional[np.ndarray]
            The B parameter of the readout filter.
        filter_gain : Optional[float]
            The gain of the readout filter.
        filter_order : Optional[int]
            The order of the readout filter.
        filter_enabled : Optional[bool]
            True if the readout filter is enabled.
        downsample_factor : Optional[int]
            Downsampling factor
        downsample_enabled : Optional[bool]
            Whether downsampler is enabled
        flux_ramp_rate_hz : float
            Flux Ramp Rate calculated from the RampMaxCnt and the digitizer
            frequency.
    """
    NUM_BANDS = 8
    CHANS_PER_BAND = 512

    def __init__(self, status):
        self.status = status

        # Reads in useful status values as attributes
        mapper_root = 'AMCc.SmurfProcessor.ChannelMapper'
        self.num_chans = self.status.get(f'{mapper_root}.NumChannels')

        # Tries to set values based on expected rogue tree
        self.mask = self.status.get(f'{mapper_root}.Mask')
        self.mask_inv = np.full((self.NUM_BANDS, self.CHANS_PER_BAND), -1)
        if self.mask is not None:
            self.mask = np.array(ast.literal_eval(self.mask))

            # Creates inverse mapping
            for i, chan in enumerate(self.mask):
                b = chan // self.CHANS_PER_BAND
                c = chan % self.CHANS_PER_BAND
                self.mask_inv[b, c] = i

        filter_root = 'AMCc.SmurfProcessor.Filter'
        self.filter_a = self.status.get(f'{filter_root}.A')
        if self.filter_a is not None:
            self.filter_a = np.array(ast.literal_eval(self.filter_a))
        self.filter_b = self.status.get(f'{filter_root}.B')
        if self.filter_b is not None:
            self.filter_b = np.array(ast.literal_eval(self.filter_b))
        self.filter_gain = self.status.get(f'{filter_root}.Gain')
        self.filter_order = self.status.get(f'{filter_root}.Order')
        self.filter_enabled = not self.status.get('{filter_root}.Disabled')

        ds_root = 'AMCc.SmurfProcessor.Downsampler'
        self.downsample_factor = self.status.get(f'{ds_root}.Factor')
        self.downsample_enabled = not self.status.get(f'{ds_root}.Disabled')

        # Tries to make resonator frequency map
        self.freq_map = np.full((self.NUM_BANDS, self.CHANS_PER_BAND), np.nan)
        band_roots = [
            f'AMCc.FpgaTopLevel.AppTop.AppCore.SysgenCryo.Base[{band}]'
            for band in range(self.NUM_BANDS)]
        for band in range(self.NUM_BANDS):
            band_root = band_roots[band]
            band_center = self.status.get(f'{band_root}.bandCenterMHz')
            subband_offset = self.status.get(f'{band_root}.toneFrequencyOffsetMHz')
            channel_offset = self.status.get(f'{band_root}.CryoChannels.centerFrequencyArray')

            # Skip band if one of these fields is None
            if None in [band_center, subband_offset, channel_offset]:
                continue

            subband_offset = np.array(ast.literal_eval(subband_offset))
            channel_offset = np.array(ast.literal_eval(channel_offset))
            self.freq_map[band] = band_center + subband_offset + channel_offset

        # Calculates flux ramp reset rate (Pulled from psmurf's code)
        rtm_root = 'AMCc.FpgaTopLevel.AppTop.AppCore.RtmCryoDet'
        ramp_max_cnt = self.status.get(f'{rtm_root}.RampMaxCnt')
        if ramp_max_cnt is None:
            self.flux_ramp_rate_hz = None
        else:
            digitizer_freq_mhz = float(self.status.get(
                f'{band_roots[0]}.digitizerFrequencyMHz', 614.4))
            ramp_max_cnt_rate_hz = 1.e6*digitizer_freq_mhz / 2.
            self.flux_ramp_rate_hz = ramp_max_cnt_rate_hz / (ramp_max_cnt + 1)

    def readout_to_smurf(self, rchan):
        """
        Converts from a readout channel number to (band, channel).

        Args
        -----
            rchans : int or List[int]
                Readout channel to convert. If a list or array is passed,
                this will return an array of bands and array of smurf channels.

        Returns
        --------
            band, channel : (int, int) or (List[int], List[int])
                The band, channel combination that is has readout channel
                ``rchan``.
        """
        abs_smurf_chan = self.mask[rchan]
        return (abs_smurf_chan // self.CHANS_PER_BAND,
                abs_smurf_chan % self.CHANS_PER_BAND)

    def smurf_to_readout(self, band, chan):
        """
        Converts from (band, channel) to a readout channel number.
        If the channel is not streaming, returns -1.

        Args:
            band : int, List[int]
                The band number, or list of band numbers corresopnding to
                channel input array.
            chan : int, List[int]
                Channel number or list of channel numbers.
        """
        return self.mask_inv[band, chan]
