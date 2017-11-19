import numpy as np


class RateCoeff():
    def get_K(self, T):
        raise NotImplementedError()


class ModifiedArrhenius(RateCoeff):
    '''
    Caculate the rate coeffcient for Modified Arrhenius Reaction
    
    ARGUMENTS:
    ==========
    b = modified arrhenius prefactor
        float
    a = arrhenius prefactor
        must be positive float
    E =  activation energy
        float
    R = ideal gas constant (optional; default 8.314)
        must be positive
     
    ATTRIBUTES:
    ===========
    b = modified arrhenius prefactor
        float
    a = arrhenius prefactor
        must be positive float
    E =  activation energy
        float
    R = ideal gas constant (optional; default 8.314)
        must be positive
    '''

    def __init__(self, a, b, E, R=8.314):
        self.b = b
        self.A = a
        self.E = E
        self.R = R

    def __repr__(self):
        return "a = {}, b = {}, E = {}, R = {}".format(self.A, self.b, self.E, self.R)

    def get_K(self, T):
        """Calculates the modified Arrhenius reaction rate coefficient
    
        INPUTS:
        =======
        A: float
           Arrhenius prefactor
           Must be positive
        b: float
           Modified Arrhenius parameter
        E: float
           Activation energy
        T: float
           Temperature
           Must be positive
        R: float, default value = 8.314
           Ideal gas constant
           Must be positive
        RETURNS:
        ========
        k: float
           Modified Arrhenius reaction rate coefficient
        """
        if self.A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if self.R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

        return self.A * T ** self.b * np.exp(-self.E / self.R / T)


class Arrhenius(RateCoeff):
    '''
    Caculate the rate coeffcient for Arrhenius Reaction
    
    ARGUMENTS:
    ==========
    a = arrhenius prefactor
        must be positive float
    E =  activation energy
        float
    R = ideal gas constant (optional; default 8.314)
        must be positive
     
    ATTRIBUTES:
    ===========
    a = arrhenius prefactor
        must be positive float
    E =  activation energy
        float
    R = ideal gas constant (optional; default 8.314)
        must be positive
    '''

    def __init__(self, a, E, R=8.314):
        self.A = a
        self.E = E
        self.R = R

    def get_K(self, T):
        """Calculates the Arrhenius reaction rate coefficient
    
        INPUTS:
        =======
        A: float
           Arrhenius prefactor
           Must be positive
        E: float
           Activation energy
        T: float
           Temperature
           Must be positive
        R: float, default value = 8.314
           Ideal gas constant
           Must be positive
        RETURNS:
        ========
        k: float
           Arrhenius reaction rate coefficient
        """

        if self.A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if self.R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

        return self.A * np.exp(-self.E / self.R / T)


class Constant(RateCoeff):
    '''
    Caculate the rate coeffcient for Modified Arrhenius Reaction
    
    ARGUMENTS:
    ==========
    k = constant rate coeffcient 
        float
        must be positive
 
     
    ATTRIBUTES:
    ===========
    k = constant rate coeffcient 
        float
        must be positive
    '''

    def __init__(self, const):
        self.k = const

    def get_K(self, T):
        """Simply returns a constant reaction rate coefficient
    
        INPUTS:
        =======
        k: float, default value = 1.0
           Constant reaction rate coefficient
        RETURNS:
        ========
        k: float
           Constant reaction rate coefficient
        """
        if self.k < 0:
            raise ValueError("Negative reaction rate coefficients are prohibited.")

        return self.k
