from pathlib import Path
from phantombuster.store import Store
from phantombuster.remoter import Scheduler
from phantombuster.remoter.serialization import register_serialisation


@register_serialisation
class Project:

    def __init__(self, path, tmp_dir=None, data_dir=None, stats_dir=None, scheduler=None):
        self.path = Path(path)
        self.scheduler = scheduler
        self._tmp_dir = tmp_dir
        self._data_dir = data_dir
        self._stats_dir = stats_dir

    def create(self):
        self.path.mkdir(exist_ok=True)

    def _get_server_path(self):
        path = self.tmp_dir / 'phantombuster.ipc'
        return f"ipc://{path}"

    def _get_server_store_path(self):
        path = self.tmp_dir / 'phantombuster.sqlite3'
        return str(path)

    def get_scheduler(self):
        if self.scheduler is None:
            self.scheduler = Scheduler(address=self._get_server_path(), start_server=True, name="scheduler", server_store_path=self._get_server_store_path())
        return self.scheduler

    @property
    def tmp_dir(self):
        if self._tmp_dir is None:
            tmp_dir = self.path / 'tmp'
        else:
            tmp_dir = Path(self._tmp_dir)
        tmp_dir.mkdir(exist_ok=True)
        return tmp_dir

    @tmp_dir.setter
    def tmp_dir(self, value):
        self._tmp_dir = value

    @property
    def data_dir(self):
        if self._data_dir is None:
            data_dir = self.path / 'data'
        else:
            data_dir = Path(self._data_dir)
        data_dir.mkdir(exist_ok=True)
        return data_dir

    @data_dir.setter
    def data_dir(self, value):
        self._data_dir = value

    @property
    def stats_dir(self):
        if self._stats_dir is None:
            stats_dir = self.path / 'stats'
        else:
            stats_dir = Path(self._stats_dir)
        stats_dir.mkdir(exist_ok=True)
        return stats_dir

    @stats_dir.setter
    def stats_dir(self, value):
        self._stats_dir = value

    @property
    def demultiplex_output_path(self):
        return str(self.data_dir / 'demultiplex.parquet')

    @property
    def demultiplex_stats_path(self):
        return str(self.stats_dir / 'demultiplex_statistics.json')

    @property
    def error_correct_output_path(self):
        return str(self.data_dir / 'error_correct.parquet')

    @property
    def error_correct_stats_path(self):
        return str(self.stats_dir / 'error_correct_statistics.json')

    @property
    def hopping_removal_output_path(self):
        return str(self.data_dir / 'hopping_removal.parquet')

    @property
    def hopping_removal_stats_path(self):
        return str(self.stats_dir / 'hopping_removal_statistics.json')

    @property
    def threshold_output_path(self):
        return str(self.data_dir / 'threshold.parquet')

    @property
    def threshold_stats_path(self):
        return str(self.stats_dir / 'threshold_statistics.json')

    def _to_json_e(self):
        return {'path': self.path, 'tmp_dir': self._tmp_dir, 'stats_dir': self._stats_dir, 'data_dir': self._data_dir}

    @classmethod
    def _from_json_e(cls, d: dict):
        return Project(**d)
