import numpy as np

from . import thermochem


class ThermochemRXNSetWrapper:
    def __init__(self, rxn_data):
        self.nasa7_coeffs = None  # todo
        nu_r, nu_p = rxn_data.get_nu()
        self.nuij = nu_p - nu_r


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

    def __init__(self, id, species, reactions, nasa):
        self.id = id
        self.reactions = reactions
        self.species = species
        self.I = len(self.species)
        self.J = len(self.reactions)
        species_set = set(self.species)
        self.species_info = {}
        ids = set()
        for s in self.species:
            self.species_info[s] = {'l': {}, 'h': {}}

            low_coeffs, low_tmin, low_tmax = nasa.get_coeffs(s, 'low')
            high_coeffs, high_tmin, high_tmax = nasa.get_coeffs(s, 'high')

            self.species_info[s]['l']['Tmax'] = low_tmax
            self.species_info[s]['l']['Tmin'] = low_tmin
            self.species_info[s]['l']['coeffs'] = low_coeffs
            self.species_info[s]['h']['Tmax'] = high_tmax
            self.species_info[s]['h']['Tmin'] = high_tmin
            self.species_info[s]['h']['coeffs'] = high_coeffs

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

        EXAMPLES:
        ========
        >>> from .parser import DataParser
        >>> from .nasa import NASACoeffs
        >>> nasa = NASACoeffs()
        >>> data_parser = DataParser()
        >>> reaction_data = data_parser.parse_file("chemkin/example_data/rxns.xml", nasa)
        >>> progress_rates = reaction_data.get_progress_rate([1,2,3,4,5,6],100)
        >>> print(progress_rates)
        [  1.06613928e-26   1.85794997e-09   1.20000000e+04]
        """
        if len(concs) != self.I:
            raise ValueError("concs must be a list of concentrations of size {}".format(self.I))
        for r in self.reactions:
            if r.reversible:
                rxnset = ThermochemRXNSetWrapper(self)
                backward = thermochem.ThermoChem(rxnset)
                kf = self.get_k(T)
                return self.__progress_rate(self.get_nu()[0], np.array(concs), backward.backward_coeffs(kf, T))
                # raise NotImplementedError("Progress rate for reversible reactions is not supported.")
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
        array: reaction rate of each specie

        EXAMPLES:
        ========
        >>> from .parser import DataParser
        >>> from .nasa import NASACoeffs
        >>> nasa = NASACoeffs()
        >>> data_parser = DataParser()
        >>> reaction_data = data_parser.parse_file("chemkin/example_data/rxns.xml", nasa)
        >>> progress_rates = reaction_data.get_progress_rate([1,2,3,4,5,6],100)
        >>> reaction_rates = reaction_data.get_reaction_rate(progress_rates)
        >>> print(reaction_rates[:4])
        [  1.20000000e+04  -1.85794997e-09  -1.20000000e+04  -1.20000000e+04]
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
