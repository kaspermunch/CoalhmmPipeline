
import os
from tables import open_file as openFile

class OpenHDF():

    def __init__(self, fileName, mode):
        self.fileName = fileName
        self.mode = mode

    def __enter__(self):
        self.hdf = openFile(self.fileName, self.mode)
        return self.hdf

    def __exit__(self, exc_type, exc_value, traceback):
        self.hdf.close()
        if os.path.exists(self.fileName) and not os.path.getsize(self.fileName):
            os.unlink(self.fileName)
