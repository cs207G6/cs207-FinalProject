import numpy as np


class RateCoeff:
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


class ReactionData:
    """
    Contains all the data related to the reaction; i.e reaction & progress rate

     ARGUMENTS:
     ==========
     id = identifier in xml
     species = reactants & product
     reactions: an array of the reactions

     ATTRIBUTES:
     ===========
     id = identifier in xml
     species = reactants & product
     reactions: an array of the reactions
    """

    def __init__(self, id, species, reactions):
        self.id = id
        self.reactions = reactions
        self.species = species
        self.I = len(self.species)
        self.J = len(self.reactions)
        species_set = set(self.species)
        ids = set()
        for r in self.reactions:
            if r.id in ids:
                raise ValueError("Duplicate id: {}".format(r.id))
            ids.add(r.id)
            for k in r.reactants:
                if k not in species_set:
                    raise ValueError("{} is not in species array.".format(k))
            for k in r.products:
                if k not in species_set:
                    raise ValueError("{} is not in species array.".format(k))

    def __len__(self):
        return self.J

    def get_nu(self):
        inv_dict = {v: k for (k, v) in enumerate(self.species)}
        nu_react = np.zeros((self.I, self.J))
        nu_prod = np.zeros((self.I, self.J))
        for (j, reaction) in enumerate(self.reactions):
            for reactant in reaction.reactants:
                nu_react[inv_dict[reactant], j] = reaction.reactants[reactant]
            for product in reaction.products:
                nu_prod[inv_dict[product], j] = reaction.products[product]
        return nu_react, nu_prod

    def get_k(self, T):
        return np.array([reaction.rate_coeff.get_K(T) for reaction in self.reactions])

    def get_progress_rate(self, concs, T):
        """Returns the progress rate of a system of irreversible, elementary reactions

        INPUTS:
        =======
        concs:    numpy array of floats
                  concentration of species
        T:        numpy array of floats
                  temperature

        RETURNS:
        ========
        omega: numpy array of floats
               size: num_reactions
               progress rate of each reaction
        """
        if len(concs) != self.I:
            raise ValueError("concs must be a list of concentrations of size {}".format(self.I))
        for r in self.reactions:
            if r.reversible:
                raise NotImplementedError("Progress rate for reversible reactions is not supported.")
            if r.type != "Elementary":
                raise NotImplementedError("Progress rate for {} reactions is not supported.", format(r.type))
        return self.__progress_rate(self.get_nu()[0], np.array(concs), self.get_k(T))

    def get_reaction_rate(self, progress_rates):
        """Returns the reaction rate of a system of irreversible, elementary reactions

        INPUTS:
        =======
        progress_rates: numpy array of floats,
                  size: num_species X num_reactions

        RETURNS:
        ========
        array: reaction rate of each species
        """
        nu_react, nu_prod = self.get_nu()
        return self.__reaction_rate(nu_react, nu_prod, progress_rates)

    def __progress_rate(self, nu_react, concs, k):
        """Returns the progress rate of a system of irreversible, elementary reactions

        INPUTS:
        =======
        nu_react: numpy array of floats,
                  size: num_species X num_reactions
                  stoichiometric coefficients for the reaction
        k:        array of floats
                  Reaction rate coefficient for the reaction
        concs:    numpy array of floats
                  concentration of species

        RETURNS:
        ========
        omega: numpy array of floats
               size: num_reactions
               progress rate of each reaction
        """
        progress = k  # Initialize progress rates with reaction rate coefficients
        for jdx, rj in enumerate(progress):
            if rj < 0:
                raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(rj))
            for idx, xi in enumerate(concs):
                nu_ij = nu_react[idx, jdx]
                if xi < 0.0:
                    raise ValueError("x{0} = {1:18.16e}:  Negative concentrations are prohibited!".format(idx, xi))
                if nu_ij < 0:
                    raise ValueError(
                        "nu_{0}{1} = {2}:  Negative stoichiometric coefficients are prohibited!".format(idx, jdx,
                                                                                                        nu_ij))

                progress[jdx] *= xi ** nu_ij
        return progress

    def __reaction_rate(self, nu_react, nu_prod, rj):
        """Returns the reaction rate of a system of irreversible, elementary reactions

        INPUTS:
        =======
        nu_react: numpy array of floats,
                  size: num_species X num_reactions
                  stoichiometric coefficients for the reactants
        nu_prod:  numpy array of floats,
                  size: num_species X num_reactions
                  stoichiometric coefficients for the products
        k:        float, default value is 10,
                  Reaction rate coefficient for the reaction
        concs:    numpy array of floats
                  concentration of species

        RETURNS:
        ========
        f: numpy array of floats
           size: num_species
           reaction rate of each specie
        """
        nu = nu_prod - nu_react
        return np.dot(nu, rj)


class Reaction:
    def __init__(self, id, reversible, type_, reactants, products, rate_coeff):
        self.id = id
        self.reversible = reversible
        self.type = type_
        self.reactants = reactants
        self.products = products
        self.rate_coeff = rate_coeff