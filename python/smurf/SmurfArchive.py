import sqlalchemy as db
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from spt3g import core
import so3g
import datetime as dt
import os
from tqdm import tqdm
import numpy as np
import yaml


num_bias_lines = 16
type_key = {
    core.G3FrameType.Observation: 0,
    core.G3FrameType.Wiring: 1,
    core.G3FrameType.Scan: 2
}

Base = declarative_base()
Session = sessionmaker()
class Frame(Base):
    """Table to store frame indexing info"""
    __tablename__ = 'frames'
    id = db.Column(db.Integer, primary_key=True)

    # location is "<path>:<index>" and this must be unique
    location = db.Column(db.String, nullable=False, unique=True)
    offset = db.Column(db.Integer, nullable=False)
    type = db.Column(db.SmallInteger, nullable=False)
    time = db.Column(db.DateTime, nullable=False)

    # Specific to data frames
    samples = db.Column(db.Integer)
    channels = db.Column(db.Integer)
    start = db.Column(db.DateTime)
    stop = db.Column(db.DateTime)

    def __repr__(self):
        return f"Frame({list(type_key.keys())[self.type]})<{self.location}>"


class Files(Base):
    """Table to store file indexing info"""
    __tablename__ = 'files'
    id = db.Column(db.Integer, primary_key=True)

    path = db.Column(db.String, nullable=False)
    start = db.Column(db.DateTime)
    stop = db.Column(db.DateTime)
    frames = db.Column(db.Integer)
    channels = db.Column(db.Integer)



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

    def add_file(self, path, session):
        """
        Indexes a single file and adds it to the sqlite database.

        Args
        ----
        path (path): Path of the file to index
        """
        reader = so3g.G3IndexedReader(path)

        total_channels = 0
        file_start, file_stop = None, None
        frame_idx = 0
        while True:
            db_frame = Frame()
            db_frame.offset = reader.Tell()

            frames = reader.Process(None)
            if not frames:
                break
            frame = frames[0]

            db_frame.location = f"{path}:{frame_idx:0>4}"
            frame_idx += 1

            if frame.type not in type_key.keys():
                continue
            db_frame.type = type_key[frame.type]
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

        db_file = Files(path=path, frames=frame_idx, channels=total_channels,
                        start=file_start, stop=file_stop)
        session.add(db_file)

    def index_archive(self, verbose=False):
        """
        Adds all files from an archive to the sqlite database.
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
                session.rollback()
                print(e)
            except RuntimeError as e:
                session.rollback()
                print(f"Failed on file {f} due to end of stream error!")
            except Exception as e:
                session.rollback()
                raise e

        session.close()

    def load_data(self, start, end, show_pb=True, load_biases=False):
        """
        Loads smurf G3 data for a given time range. For the specified time range
        this will return a chunk of data that includes that time range.

        Args
        -----
        start (timestamp): start timestamp
        end   (timestamp): end timestamp
        show_pb (bool, optional): If True, will show progress bar.
        """
        session = self.Session()
        frames = session.query(Frame).filter(
            Frame.type == type_key[core.G3FrameType.Scan],
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

        timestamps = np.full((samples,), np.nan)
        data = np.full((channels, samples), 0, dtype=np.int32)
        if load_biases:
            biases = np.full((num_bias_lines, samples), 0, dtype=np.int32)
        else:
            baises = None

        cur_sample = 0
        cur_file = None
        for frame_info in tqdm(frames, total=num_frames, disable=(not show_pb)):
            file = frame_info.location.split(':')[0]
            if file != cur_file:
                reader = so3g.G3IndexedReader(file)
                cur_file = file

            reader.Seek(frame_info.offset)
            frame = reader.Process(None)[0]
            nsamp = frame['data'].n_samples
            timestamps[cur_sample:cur_sample + nsamp] = frame['data'].times()
            key_order = [int(k[1:]) for k in frame['data'].keys()]
            data[key_order, cur_sample:cur_sample + nsamp] = frame['data']
            if load_biases:
                bias_order = [int(k[-2:]) for k in frame['tes_biases'].keys()]
                biases[bias_order, cur_sample:cur_sample + nsamp] = frame['tes_biases']
            cur_sample += nsamp

        timestamps /= core.G3Units.s
        return timestamps, data, biases

    def load_status(self, time):
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
            Frame.type == type_key[core.G3FrameType.Observation],
            Frame.time <= dt.datetime.fromtimestamp(time)
        ).order_by(Frame.time.desc()).first()

        status_frames = session.query(Frame).filter(
            Frame.type == type_key[core.G3FrameType.Wiring],
            Frame.time >= session_start,
            Frame.time <= dt.datetime.fromtimestamp(time)
        ).order_by(Frame.time)

        status = {}
        cur_file = None
        for frame_info in tqdm(status_frames.all()):
            file = frame_info.location.split(':')[0]
            if file != cur_file:
                reader = so3g.G3IndexedReader(file)
                cur_file = file
            reader.Seek(frame_info.offset)
            frame = reader.Process(None)[0]
            status.update(yaml.safe_load(frame['status']))

        return status

