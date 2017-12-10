import base64

import h5py
import matplotlib.pyplot as plt
from io import BytesIO


class TimeEvo:
    def __init__(self, file):
        self.file = h5py.File(file, 'r')
        self.scenarios = list(self.file.keys())

    def plot(self, scenario):
        data = self.file[scenario + '/truth']
        time = self.file[scenario + '/time']
        pic_width = 1200 // 75
        pic_length = 800 // 75
        plt.figure(figsize=(pic_width, pic_length))
        plt.plot(time.value, data.value[:, -1], label="Temperature")
        plt.title("Temperature Evolution")
        plt.legend()
        figfile = BytesIO()
        plt.savefig(figfile, format='png')
        figfile.seek(0)  # rewind to beginning of file
        figdata_png = base64.b64encode(figfile.getvalue())
        return figdata_png.decode('utf8')
