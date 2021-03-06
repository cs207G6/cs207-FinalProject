import numpy as np


class RateCoeff:
    def get_K(self, T):
        raise NotImplementedError()


class ModifiedArrhenius(RateCoeff):
    """
    Calculate the rate coefficient for Modified Arrhenius Reaction

    Attributes
    ----------
    b: float
     modified arrhenius prefactor
    a: float
     arrhenius prefactor
    E: float
      activation energy
    R: float
     ideal gas constant (optional; default 8.314). must be positive
    """

    def __init__(self, a, b, E, R=8.314):
        """
        Create a new instance of ModifiedArrhenius

        Parameters
        ----------
        b: float
         modified arrhenius prefactor
        a: float
         arrhenius prefactor,must be positive float
        E: float
          activation energy
        R: float
         ideal gas constant (optional; default 8.314). must be positive
        """
        self.b = b
        self.A = a
        self.E = E
        self.R = R

    def get_K(self, T):
        """
        Calculates the modified Arrhenius reaction rate coefficient
    
        Parameters
        ----------
        T: float
           Temperature
           Must be positive

        Returns
        -------
        k: float
           Modified Arrhenius reaction rate coefficient

        Examples
        --------
        >>> ModifiedArrhenius(10,20,30).get_K(50)
        8.872748484824047e+34
        """
        if self.A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if self.R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

        return self.A * T ** self.b * np.exp(-self.E / self.R / T)


class Arrhenius(RateCoeff):
    """
    Calculate the rate coefficient for Arrhenius Reaction

    Attributes
    ----------
    a: float
     arrhenius prefactor. must be positive float
    E: float
      activation energy
    R: float
     ideal gas constant (optional; default 8.314). must be positive
    """

    def __init__(self, a, E, R=8.314):
        """
        Create a new Arrhenius instance

        Parameters
        ----------
        a: float
         arrhenius prefactor. must be positive float
        E: float
         activation energy
        R: float
         ideal gas constant (optional; default 8.314). must be positive
        """
        self.A = a
        self.E = E
        self.R = R

    def get_K(self, T):
        """
        Calculates the Arrhenius reaction rate coefficient
    
        Parameters
        ----------
        T: float
           Temperature
           Must be positive

        Returns
        -------
        k: float
           Arrhenius reaction rate coefficient

        Examples
        --------
        >>> Arrhenius(10,20,30).get_K(50)
        9.8675516180719569
        """

        if self.A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if self.R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

        return self.A * np.exp(-self.E / self.R / T)


class Constant(RateCoeff):
    """
    Get constant rate coefficient

    Attributes
    ----------
    k: float
        constant rate coefficient, must be positive
    """

    def __init__(self, const):
        """
        Create a new instance for constant rate coefficient

        Parameters
        ----------
        const: float
            constant rate coefficient, must be positive
        """
        self.k = const

    def get_K(self, T):
        """
        Simply returns a constant reaction rate coefficient

        Parameters
        ----------
        T: float
           Temperature
           Must be positive

        Returns
        -------
        k: float
           constant reaction rate coefficient

        Examples
        --------
        >>> Constant(1).get_K(50)
        1
        """
        if self.k < 0:
            raise ValueError("Negative reaction rate coefficients are prohibited.")

        return self.k
