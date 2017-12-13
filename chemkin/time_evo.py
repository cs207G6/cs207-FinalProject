import base64

import h5py
import matplotlib
from io import BytesIO

matplotlib.use('Agg')


class TimeEvo:
    def __init__(self, file):
        """
        Create a TimeEvo instance for plotting temperature evolution

        Parameters
        ----------
        file: str
            path to the hdf5 file
        """
        self.file = h5py.File(file, 'r')
        self.scenarios = list(self.file.keys())

    def plot(self, scenario, pic_width=16, pic_length=10):
        """
        Returns a base64 encoded image of time temperature evolution of specified scenario

        Parameters
        ----------
        scenario: str
            name of the scenario
        pic_width:
            width of output picture
        pic_length:
            height of output picture

        Returns
        -------
        str
            base64 encoded image (png)

        Examples
        --------
        >>> time_evo = TimeEvo("chemkin/example_data/detailed_profile.h5")
        >>> plot = time_evo.plot(time_evo.scenarios[0], 0.1, 0.1)
        """
        import matplotlib.pyplot as plt
        data = self.file[scenario + '/truth']
        time = self.file[scenario + '/time']
        plt.figure(figsize=(pic_width, pic_length))
        plt.plot(time.value, data.value[:, -1], label="Temperature")
        plt.xlabel("Time")
        plt.ylabel("Temp")
        plt.title("Temperature Evolution")
        plt.legend()
        figfile = BytesIO()
        plt.savefig(figfile, format='png')
        figfile.seek(0)  # rewind to beginning of file
        figdata_png = base64.b64encode(figfile.getvalue())
        return figdata_png.decode('utf8')
