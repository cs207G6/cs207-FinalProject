import base64

import h5py
import matplotlib.pyplot as plt
from io import BytesIO


class TimeEvo:
    def __init__(self, file):
        self.file = h5py.File(file, 'r')
        self.scenarios = list(self.file.keys())

    def plot(self, scenario, pic_width=16, pic_length=10):
        """

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
        >>> time_evo.plot(time_evo.scenarios[0], 0.1, 0.1)
        'iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAPYQAAD2EBqD+naQAAANJJREFUGJWFkDFOw0AQRZ+XtSxZ2m1oseg4AU0OgMQJ0tBt4dNQpqHmGG4pKGjjCwxyk8q7JFIIm6EIcUEiGGmKP3oavZki56zDMOCcoygKfpeqklICEVHg37bOOQBEhPcPZf70CsD89gqA7WbNor3DppQA8N5T5i9MVXN9WfP4MAMgxsiiBdM0zeSz3e0BqKw5cTUiMoXPfADLizOg934Ke9XD8Mz1NsYIQPv8xnK1+wFPOOzR8WUpmKqmKg33Nw3HBeM4EkKAvu///F8IQbuu028NBl10rRrWqQAAAABJRU5ErkJggg=='
        """
        data = self.file[scenario + '/truth']
        time = self.file[scenario + '/time']
        plt.figure(figsize=(pic_width, pic_length))
        plt.plot(time.value, data.value[:, -1], label="Temperature")
        plt.title("Temperature Evolution")
        plt.legend()
        figfile = BytesIO()
        plt.savefig(figfile, format='png')
        figfile.seek(0)  # rewind to beginning of file
        figdata_png = base64.b64encode(figfile.getvalue())
        return figdata_png.decode('utf8')
