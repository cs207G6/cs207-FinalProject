import numpy as np
class ReactionData():
    '''
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
    '''
    def __init__(self, id, species, reactions):
        self.id = id
        self.reactions = reactions
        self.species = species
        self.I = len(self.species)
        self.J = len(self.reactions)
        species_set = set(self.species)
        for r in self.reactions:
            for k in r.reactants:
                if k not in species_set:
                    raise ValueError("{} is not in species array.".format(k))
            for k in r.products:
                if k not in species_set:
                    raise ValueError("{} is not in species array.".format(k))
        
    def get_nu(self):
        inv_dict = {v:k for (k,v) in enumerate(self.species)}
        nu_react = np.zeros((self.I, self.J))
        nu_prod = np.zeros((self.I, self.J))
        for (j,reaction) in enumerate(self.reactions):
            for reactant in reaction.reactants:
                nu_react[inv_dict[reactant],j]=reaction.reactants[reactant]
            for product in reaction.products:
                nu_prod[inv_dict[product],j]=reaction.products[product]
        return nu_react, nu_prod
    
    def get_k(self, T):
        return np.array([reaction.rate_coeff.get_K(T) for reaction in self.reactions])
    
    def get_progress_rate(self, concs, T):
        if len(concs) != self.I:
            raise ValueError("concs must be a list of concentrations of size {}".format(self.I))
        for r in self.reactions:
            if r.reversible:
                raise NotImplementedError("Progress rate for reversible reactions is not supported.")
            if r.type != "Elementary":
                raise NotImplementedError("Progress rate for {} reactions is not supported.",format(r.type))
        return self.__progress_rate(self.get_nu()[0], np.array(concs), self.get_k(T))
    
    def get_reaction_rate(self, reaction_rates):
        nu_react, nu_prod = self.get_nu()
        return self.__reaction_rate(nu_react, nu_prod, reaction_rates)
    
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

        EXAMPLES:
        =========
        >>> progress_rate_2(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([2.0, 1.0, 1.0]), 10.0)
        array([ 40.,  20.])
        """
        progress = k # Initialize progress rates with reaction rate coefficients
        for jdx, rj in enumerate(progress):
            if rj < 0:
                raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(rj))
            for idx, xi in enumerate(concs):
                nu_ij = nu_react[idx,jdx]
                if xi  < 0.0:
                    raise ValueError("x{0} = {1:18.16e}:  Negative concentrations are prohibited!".format(idx, xi))
                if nu_ij < 0:
                    raise ValueError("nu_{0}{1} = {2}:  Negative stoichiometric coefficients are prohibited!".format(idx, jdx, nu_ij))

                progress[jdx] *= xi**nu_ij
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

        EXAMPLES:
        =========
        >>> nu_react = np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
        >>> nu_prod = np.array([[0.0, 1.0], [0.0, 2.0], [1.0, 0.0]])
        >>> r = np.array([ 40.,  20.])
        >>> reaction_rate(nu_react, nu_prod, r)
        array([-20., -40.,   0.])
        """
        nu = nu_prod - nu_react
        return np.dot(nu, rj)
        
class RateCoeff():
    def get_K(self, T):
        raise NotImplementedError()

class ModifiedArrhenius(RateCoeff):
    def __init__(self, a, b, E, R = 8.314):
        self.b = b
        self.A = a
        self.E = E
        self.R = R
    
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

        EXAMPLES:
        =========
        >>> k_mod_arr(2.0, -0.5, 3.0, 100.0)
        0.19927962618542916
        """
        if self.A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if self.R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

        return self.A * T**self.b * np.exp(-self.E / self.R / T)
    
class Arrhenius(RateCoeff):
    def __init__(self, a, E,R=8.314):
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

        EXAMPLES:
        =========
        >>> k_arr(2.0, 3.0, 100.0)
        1.9927962618542914
        """

        if self.A < 0.0:
            raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

        if T < 0.0:
            raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(T))

        if self.R < 0.0:
            raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

        return self.A * np.exp(-self.E / self.R / T)
    
class Constant(RateCoeff):
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

        EXAMPLES:
        =========
        >>> k_const(5.0)
        5.0
        """
        if self.k < 0:
            raise ValueError("Negative reaction rate coefficients are prohibited.")

        return self.k 
        
class Reaction():
    def __init__(self, id, reversible, type_, reactants, products, rate_coeff):
        self.id = id
        self.reversible = reversible
        self.type = type_
        self.reactants = reactants
        self.products = products
        self.rate_coeff = rate_coeff

class DataParser():
    def __init__(self):
        pass
    
    def _parse_rate_coeff(self, root):
        root = root.getchildren()[0]
        if root.tag == "Arrhenius":
            A = float(root.find('A').text)
            E = float(root.find('E').text)
            if root.find('R') is not None:
                return Arrhenius(a=A, E=E, R=float(R.text))
            return Arrhenius(a=A, E=E)
        elif root.tag == "modifiedArrhenius":
            A = float(root.find('A').text)
            b = float(root.find('b').text)
            E = float(root.find('E').text)
            if root.find('R') is not None:
                return ModifiedArrhenius(a=A, b=b, E=E, R=float(R.text))
            return ModifiedArrhenius(a=A, b=b, E=E)
        elif root.tag == "Constant":
            k = float(root.find('k').text)
            return Constant(k)
        else:
            raise NotImplementedError("{} rate coefficient is not supported".format(root.tag))
    
    def _get_bool(self, v):
        v = str(v).lower()
        if v == 'yes' or v == 'true' or v == '1':
            return True
        else:
            return False
    
    def _get_reactants(self, root):
        s = root.text.strip()
        reactants = s.split()
        result = {}
        for r in reactants:
            parts = r.split(':')
            result[parts[0]]=parts[1]
        return result
    
    def _parse_reaction(self, root):
        rate_coeff = self._parse_rate_coeff(root.find('rateCoeff'))
        reversible = self._get_bool(root.get('reversible'))
        type_ = root.get('type')
        id = root.get('id')
        reactants = self._get_reactants(root.find('reactants'))
        products = self._get_reactants(root.find('products'))
        return Reaction(id, reversible=reversible, type_=type_, reactants=reactants, products=products, rate_coeff=rate_coeff)
    
    def parse_file(self, filename):
        import xml.etree.ElementTree as ET
        tree = ET.parse(filename)
        id = tree.find("reactionData").get("id")
        species = []
        for child in tree.find('phase'):
            species += child.text.split()
        reactions = []
        for (i,r) in enumerate(tree.find('reactionData').findall('reaction')):
            reaction = self._parse_reaction(r)
            reactions.append(reaction)
        return ReactionData(id, species, reactions)

