import numpy as np

from . import thermochem


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
        self.nasa = {}
        ids = set()
        for s in self.species:
            self.nasa[s] = {'l': {}, 'h': {}}

            low_coeffs, low_tmin, low_tmax = nasa.get_coeffs(s, 'low')
            high_coeffs, high_tmin, high_tmax = nasa.get_coeffs(s, 'high')

            self.nasa[s]['l']['Tmax'] = low_tmax
            self.nasa[s]['l']['Tmin'] = low_tmin
            self.nasa[s]['l']['coeffs'] = low_coeffs
            self.nasa[s]['h']['Tmax'] = high_tmax
            self.nasa[s]['h']['Tmin'] = high_tmin
            self.nasa[s]['h']['coeffs'] = high_coeffs

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

    def get_nasa_coeff(self, species, temp):
        if species not in self.nasa:
            raise KeyError("NASA coefficient for {} is not specified".format(species))
        nasa = self.nasa[species]
        if nasa['l']['Tmin'] <= temp <= nasa['l']['Tmax']:
            nasa_coeff = nasa['l']['coeffs']
        elif nasa['h']['Tmin'] <= temp <= nasa['h']['Tmax']:
            nasa_coeff = nasa['h']['coeffs']
        else:
            raise ValueError("NASA coefficient for {} at T={} is not specified".format(species, temp))
        return np.array(nasa_coeff)

    def get_nasa_coeff_matrix(self, T):
        inv_dict = {v: k for (k, v) in enumerate(self.species)}
        result = np.zeros((self.I, 7))
        species_set = set()

        for reaction in self.reactions:
            if reaction.reversible:
                species_set = species_set | set(reaction.reactants.keys())
                species_set = species_set | set(reaction.products.keys())
        for species in species_set:
            result[inv_dict[species], :] = self.get_nasa_coeff(species, T)

        return result

    def get_k(self, T):
        return np.array([reaction.rate_coeff.get_K(T) for reaction in self.reactions])

    def get_kb(self, kf, nu, T):
        nasa = self.get_nasa_coeff_matrix(T)
        tc = thermochem.ThermoChem(thermochem.ThermochemRXNSetWrapper(nasa), T)
        result = [0] * len(self.reactions)
        for j, reaction in enumerate(self.reactions):
            if reaction.reversible:
                result[j] = tc.backward_coeffs(nu[:, j], kf[j])
            else:
                pass
        return np.array(result)

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
            if r.type != "Elementary":
                raise NotImplementedError("Progress rate for {} reactions is not supported.", format(r.type))
        nus = self.get_nu()

        nu = nus[1] - nus[0]

        kf = self.get_k(T)
        forward_part = self.__progress_rate(self.get_nu()[0], np.array(concs), kf)

        kb = self.get_kb(kf=kf, nu=nu, T=T)
        backward_part = self.__progress_rate(self.get_nu()[1], np.array(concs), kb)

        return forward_part - backward_part

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
