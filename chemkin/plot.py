# coding: utf-8

# In[1]:

import base64
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np


def range_data_collection(user_data, input_concentration, lower_T, upper_T, current_T):
    """Returns all processed data required for the progress/reaction rates plots

        INPUTS:
        =======
        user_data:              ReactionData class object
                                data of all elementary reactions from website that user uploaded
        input_concentration:    python list of floats or integers
                                The concentrations of each species given by users
        lower_T:                float or integer
                                The lower bound of temperature range
        upper_T:                float or integer
                                the upper bound of temperature range
        current_T:              float or integer
                                The temperature of the current reaction
        
        RETURNS:
        ========
        temp_range:             python list of 100 floats
                                100 temperatures evenly spaced within the given temperature range
        progress_rates_list:    python list of  100 numpy arrays
                                100 progress rate arrays for 100 temperatures to be plotted respectively. Each array contains n floats 
                                where n is the number of elementary reactions in the given system
        reaction_rates_list:    python list of  100 numpy arrays
                                100 reaction rate arrays for 100 temperatures to be plotted respectively. Each array contains n floats 
                                where n is the number of species in the given system
        current_T:              float or integer
                                The temperature of the current reaction
        species:                python list of strings
                                names of all species in the given system
        progress_rates_current: numpy array
                                the progress rates for all elementary reactions of at the given current temperature
        reaction_rates_current: numpy array
                                the reaction rates for all species of at the given current temperature

        EXAMPLES:
        ========
        >>> from .parser import DataParser
        >>> from .nasa import NASACoeffs
        >>> nasa = NASACoeffs() 
        >>> data_parser = DataParser() 
        >>> user_data = data_parser.parse_file('chemkin/example_data/rxns.xml', nasa)
        >>> input_concentration = [10,30,30,40,50,60]
        >>> lower_T = 1000
        >>> upper_T = 2000
        >>> current_T = 150
        >>> pic_width = 0.1
        >>> pic_length = 0.1
        >>> (temp_range, progress_rates_list, reaction_rates_list, current_T, species, progress_rates_current, reaction_rates_current)=range_data_collection(user_data, input_concentration, lower_T, upper_T, current_T)
    """

    reaction_data = user_data
    species = reaction_data.species
    concentration = input_concentration
    temp_range = list(np.linspace(lower_T, upper_T, num=100))
    progress_rates_list = [None] * len(temp_range)
    reaction_rates_list = [None] * len(temp_range)

    for i, temp in enumerate(temp_range):
        progress_rates = reaction_data.get_progress_rate(concentration, temp)
        progress_rates_list[i] = progress_rates
        reaction_rates_list[i] = reaction_data.get_reaction_rate(progress_rates)
    # print(len(temp_range),len(progress_rates_list))

    progress_rates_current = reaction_data.get_progress_rate(concentration, current_T)
    reaction_rates_current = reaction_data.get_reaction_rate(progress_rates_current)
    return temp_range, progress_rates_list, reaction_rates_list, current_T, species, progress_rates_current, reaction_rates_current


def progress_rate_plot_generation(T_range, progress_rate_range, current_T, progress_rates_current, pic_width,
                                  pic_length):
    """Returns the base64-encoded plot for progress rates of all elementary reacitons vs temperature range

        INPUTS:
        =======
        T_range:                python list of 100 floats
                                100 temperatures evenly spaced within the given temperature range
        progress_rates_range:   python list of  100 numpy arrays
                                100 progress rate arrays for 100 temperatures to be plotted respectively. Each array contains n floats 
                                where n is the number of elementary reactions in the given system
        current_T:              float or integer
                                The temperature of the current reaction
        progress_rates_current: numpy array
                                the progress rates for all elementary reactions of at the given current temperature
        pic_width:              float
                                the desired width of the output plot image(pixel/dpi)
        pic_length:             float
                                the desired lengthth of the output plot image(pixel/dpi)                       
        
        RETURNS:
        ========
        png_str:                string
                                the base64-encoded png image of the progress rates vs temperature plot

        EXAMPLES:
        ========
        >>> from .parser import DataParser
        >>> from .nasa import NASACoeffs
        >>> nasa = NASACoeffs() 
        >>> data_parser = DataParser() 
        >>> user_data = data_parser.parse_file('chemkin/example_data/rxns.xml', nasa)
        >>> input_concentration = [10,30,30,40,50,60]
        >>> lower_T = 1000
        >>> upper_T = 2000
        >>> current_T = 150
        >>> pic_width = 0.1
        >>> pic_length = 0.1
        >>> (temp_range, progress_rates_list, reaction_rates_list, current_T, species, progress_rates_current, reaction_rates_current)=range_data_collection(user_data, input_concentration, lower_T, upper_T, current_T)
        >>> img = progress_rate_plot_generation(temp_range, progress_rates_list, current_T, progress_rates_current, pic_width,pic_length)
    """

    x = T_range
    curr_T = current_T

    # generate colors for the scatterplot
    # currT_index = x.index(curr_T)
    # colors = ['blue' for i in range(currT_index)]
    # colors.append('red')
    # colors.extend(['blue' for i in range(currT_index+1,len(x))])


    # generate plot
    plt.figure(figsize=(pic_width, pic_length))

    reaction_num = len(progress_rate_range[0])
    for i in range(reaction_num):
        # generate a curve for each elementary reaction
        y = [e[i] for e in progress_rate_range]
        plt.plot(x, y,alpha=0.6)
        plt.scatter(x, y, label='Elementary Reaction ' + str(i + 1),alpha=0.6)
        if i == 0:
            plt.plot(curr_T, progress_rates_current[i], '^r',label='Current Temperature',markersize=12)
        else:
            plt.plot(curr_T, progress_rates_current[i], '^r',markersize=12)
    plt.xlabel("Temperature")
    plt.ylabel("Progress Rate")
    plt.title("Progress Rate vs Temperature by Reactions")
    plt.legend()

    # output plot in base64 format
    figfile = BytesIO()
    plt.savefig(figfile, format='png')
    figfile.seek(0)  # rewind to beginning of file
    figdata_png = base64.b64encode(figfile.getvalue())
    png_str = figdata_png.decode('utf8')
    return png_str


def reaction_rate_plot_generation(T_range, reaction_rate_range, current_T, reaction_rates_current, species, pic_width,
                                  pic_length):
    """Returns the base64-encoded plot for reaction rates of all species vs temperature range

        INPUTS:
        =======
        T_range:                python list of 100 floats
                                100 temperatures evenly spaced within the given temperature range
        reaction_rates_range:   python list of  100 numpy arrays
                                100 reaction rate arrays for 100 temperatures to be plotted respectively. Each array contains n floats 
                                where n is the number of species in the given system
        current_T:              float or integer
                                The temperature of the current reaction
        reaction_rates_current: numpy array
                                the reaction rates for all species of at the given current temperature
        species:                python list of strings
                                names of all species in the given system
        pic_width:              float
                                the desired width of the output plot image(pixel/dpi)
        pic_length:             float
                                the desired length of the output plot image(pixel/dpi)                     
        
        RETURNS:
        ========
        png_str:                string
                                the base64-encoded png image of the reaction rates vs temperature plot

        EXAMPLES:
        ========
        >>> from .parser import DataParser
        >>> from .nasa import NASACoeffs
        >>> nasa = NASACoeffs() 
        >>> data_parser = DataParser() 
        >>> user_data = data_parser.parse_file('chemkin/example_data/rxns.xml', nasa)
        >>> input_concentration = [10,30,30,40,50,60]
        >>> lower_T = 1000
        >>> upper_T = 2000
        >>> current_T = 150
        >>> pic_width = 0.1
        >>> pic_length = 0.1
        >>> (temp_range, progress_rates_list, reaction_rates_list, current_T, species, progress_rates_current, reaction_rates_current)=range_data_collection(user_data, input_concentration, lower_T, upper_T, current_T)
        >>> img = reaction_rate_plot_generation(temp_range, reaction_rates_list, current_T, reaction_rates_current, species, pic_width,pic_length)
    """
    x = T_range
    curr_T = current_T

    # generate colors for the scatterplot
    # currT_index = x.index(curr_T)
    # colors = ['blue' for i in range(currT_index)]
    # colors.append('red')
    # colors.extend(['blue' for i in range(currT_index+1,len(x))])


    # generate plot
    plt.figure(figsize=(pic_width, pic_length))

    species_num = len(reaction_rate_range[0])
    for i in range(species_num):
        # generate a curve for each elementary reaction
        y = [e[i] for e in reaction_rate_range]
        plt.plot(x, y, label=None,alpha=0.6)
        plt.scatter(x, y, label=species[i],alpha=0.6)
        if i == 0:
            plt.plot(curr_T, reaction_rates_current[i], '^r',label='Current Temperature',markersize=12)
        else:
            plt.plot(curr_T, reaction_rates_current[i], '^r',markersize=12)
    plt.xlabel("Temperature")
    plt.ylabel("Reaction Rate")
    plt.title("Reaction Rate vs Temperature by Species")
    plt.legend()

    # output plot in base64 format
    figfile = BytesIO()
    plt.savefig(figfile, format='png')
    figfile.seek(0)  # rewind to beginning of file
    figdata_png = base64.b64encode(figfile.getvalue())
    png_str = figdata_png.decode('utf8')
    return png_str


