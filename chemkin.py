import numpy as np
import sqlite3
import numpy as np
import pandas as pd

class ThermochemRXNSetWrapper:
    def __init__(self, rxnData:ReactionData):
        self.nasa7_coeffs = None    # todo
        nu_r, nu_p = rxnData.get_nu()
        self.nuij = nu_p - nu_r

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
        >>> data_parser = DataParser()
        >>> reaction_data = data_parser.parse_file("data/rxns.xml")
        >>> progress_rates = reaction_data.get_progress_rate([1,2,3,4,5,6],100)
        >>> print(progress_rates)
        [  1.06613928e-26   1.85794997e-09   1.20000000e+04]
        """
        if len(concs) != self.I:
            raise ValueError("concs must be a list of concentrations of size {}".format(self.I))
        for r in self.reactions:
            #if r.reversible:
                #backward = thermochem(rxnset) #???what is rxnset
                #return self.__progress_rate(self.get_nu()[0], np.array(concs), backward.backward_coeffs(kf,T)) # ???what is kf?
                raise NotImplementedError("Progress rate for reversible reactions is not supported.")
            if r.type != "Elementary":
                raise NotImplementedError("Progress rate for {} reactions is not supported.",format(r.type))
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
        >>> data_parser = DataParser()
        >>> reaction_data = data_parser.parse_file("data/rxns.xml")
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
        """
        nu = nu_prod - nu_react
        return np.dot(nu, rj)
        
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
    def __init__(self, a, b, E, R = 8.314):
        self.b = b
        self.A = a
        self.E = E
        self.R = R
    
    def __repr__(self):
         return "a = {}, b = {}, E = {}, R = {}".format(self.A, self.b,self.E,self.R)
    
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

        return self.A * T**self.b * np.exp(-self.E / self.R / T)
    
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
        
class Reaction():
    def __init__(self, id, reversible, type_, reactants, products, rate_coeff):
        self.id = id
        self.reversible = reversible
        self.type = type_
        self.reactants = reactants
        self.products = products
        self.rate_coeff = rate_coeff

class thermochem:
    """Methods for calculating the backward reaction rate.

    Cp_over_R: Returns specific heat of each specie given by 
               the NASA polynomials.
    H_over_RT:  Returns the enthalpy of each specie given by 
                the NASA polynomials.
    S_over_R: Returns the entropy of each specie given by 
              the NASA polynomials.
    backward_coeffs:  Returns the backward reaction rate 
                      coefficient for reach reaction.

    Please see the notes in each routine for clarifications and 
    warnings.  You will need to customize these methods (and 
    likely the entire class) to suit your own code base.  
    Nevertheless, it is hoped that you will find these methods 
    to be of some use.
    """

    def __init__(self, rxnset):
        self.rxnset = rxnset
        self.p0 = 1.0e+05 # Pa
        self.R = 8.3144598 # J / mol / K
        self.gamma = np.sum(self.rxnset.nuij, axis=0)

    def Cp_over_R(self, T):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.nasa7_coeffs

        Cp_R = (a[:,0] + a[:,1] * T + a[:,2] * T**2.0 
                + a[:,3] * T**3.0 + a[:,4] * T**4.0)

        return Cp_R

    def H_over_RT(self, T):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.nasa7_coeffs

        H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0 
                + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0 
                + a[:,5] / T)

        return H_RT
               

    def S_over_R(self, T):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.nasa7_coeffs

        S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0 
               + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])

        return S_R

    def backward_coeffs(self, kf, T):

        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = np.dot(self.rxnset.nuij.T, self.H_over_RT(T))
        delta_S_over_R = np.dot(self.rxnset.nuij.T, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction 
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = self.p0 / self.R / T

        # Ke
        kb = fact**self.gamma * np.exp(delta_G_over_RT)

        return kf / kb



class DataParser():
    """
    XML parser for reaction data xml files.
    """
    
    def __init__(self):
        pass
    
    def _parse_rate_coeff(self, root):
        """
        Parse 'rateCoeff' element and return RateCoeff 
    
        Arguments
        ----------
        root: Element
            rateCoeff Element 
            
        Return
        ----------
        RateCoeff
            RateCoeff object
        """
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
        """
        Convert an xml string boolean to python bool
        
        Arguments
        ----------
        v: str
            xml string boolean
        
        Return
        ----------
        bool
            boolean value
        """
        v = str(v).lower()
        if v == 'yes' or v == 'true' or v == '1':
            return True
        else:
            return False
    
    def _get_reactants(self, root):
        """
        Parse 'reactants' or 'products' element and return dictionary of reactants/products 
    
        Arguments
        ----------
        root: Element
            reactants or products Element 
            
        Return
        ----------
        dict
            a string to float array where keys are name of reactants/products and value being their Stoichiometric coefficients
        """
        s = root.text.strip()
        reactants = s.split()
        result = {}
        for r in reactants:
            parts = r.split(':')
            result[parts[0]]=float(parts[1])
        return result
    
    def _parse_reaction(self, root):
        """
        Parse 'reaction' element and return Reaction 
    
        Arguments
        ----------
        root: Element
            reaction Element 
            
        Return
        ----------
        Reaction
            Reaction object
        """
        rate_coeff = self._parse_rate_coeff(root.find('rateCoeff'))
        reversible = self._get_bool(root.get('reversible'))
        type_ = root.get('type')
        id = root.get('id')
        reactants = self._get_reactants(root.find('reactants'))
        products = self._get_reactants(root.find('products'))
        return Reaction(id, reversible=reversible, type_=type_, reactants=reactants, products=products, rate_coeff=rate_coeff)
    
    def parse_file(self, filename):
        """
        Parse a reaction xml file and return ReactionData object
    
        Arguments
        ----------
        filename: str
            filename of reaction xml file 
            
        Return
        ----------
        ReactionData
            parsed ReactionData object
        """
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
    
class NASACoeffs(self,species):   
     def create_db(self,filename):
            """
            Parse a reaction xml file and return sql database with 

            Arguments
            ----------
            filename: str
                filename of reaction xml file 

            Return
            ----------
            HW10_demo.sqlite
                database with Nasa coefficient for each species
            """
            #Import SQLlite
            pd.set_option('display.width', 500)
            pd.set_option('display.max_columns', 100)
            pd.set_option('display.notebook_repr_html', True)

            db = sqlite3.connect('HW10_demo.sqlite')
            cursor = db.cursor()
            cursor.execute("DROP TABLE IF EXISTS LOW")
            cursor.execute("DROP TABLE IF EXISTS HIGH")
            cursor.execute("PRAGMA foreign_keys=1")

            #Create High and Low tables
            cursor.execute('''CREATE TABLE LOW ( 
                           SPECIES_NAME TEXT NOT NULL, 
                           TLOW REAL NOT NULL, 
                           THIGH REAL NOT NULL, 
                           COEFF_1 TEXT NOT NULL,
                           COEFF_2 TEXT NOT NULL,
                           COEFF_3 TEXT NOT NULL,
                           COEFF_4 TEXT NOT NULL,
                           COEFF_5 TEXT NOT NULL,
                           COEFF_6 TEXT NOT NULL,
                           COEFF_7 TEXT NOT NULL)''')

            # Commit changes to the database
            db.commit()
            cursor.execute('''CREATE TABLE HIGH (
                           SPECIES_NAME TEXT NOT NULL, 
                           TLOW REAL NOT NULL, 
                           THIGH REAL NOT NULL, 
                           COEFF_1 TEXT NOT NULL,
                           COEFF_2 TEXT NOT NULL,
                           COEFF_3 TEXT NOT NULL,
                           COEFF_4 TEXT NOT NULL,
                           COEFF_5 TEXT NOT NULL,
                           COEFF_6 TEXT NOT NULL,
                           COEFF_7 TEXT NOT NULL)''')
            db.commit()
            
            #Parse XML to get info for each species
            tree = ET.parse(filename)
            root = tree.getroot()

            #get species
            species = root.find('speciesData').findall('species')

            for specie in species:
                name = specie.get('name')

                #get low temp high/low and coeffs for each specie
                NASA = specie.find('thermo').findall('NASA')

                #get low info
                low_tmax = NASA[0].get('Tmax')
                low_tmin = NASA[0].get('Tmin')
                Low_C_1,Low_C_2,Low_C_3,Low_C_4,Low_C_5,Low_C_6,Low_C_7 = NASA[0].find('floatArray').text.split()
                lows_to_insert = (name,float(low_tmin),float(low_tmax),Low_C_1,Low_C_2,Low_C_3,Low_C_4,Low_C_5,Low_C_6,Low_C_7)


                #get low info
                high_tmax = NASA[1].get('Tmax')
                high_tmin = NASA[1].get('Tmin')
                High_C_1,High_C_2,High_C_3,High_C_4,High_C_5,High_C_6,High_C_7 = NASA[1].find('floatArray').text.split()
                high_to_insert = name,float(high_tmin),float(high_tmax),High_C_1,High_C_2,High_C_3,High_C_4,High_C_5,High_C_6,High_C_7

                #Insert the values for each species into table
                cursor.execute('''INSERT INTO LOW 
                              (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', lows_to_insert)
                cursor.execute('''INSERT INTO HIGH 
                              (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7)
                              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', high_to_insert)

    def get_coeffs(species_name, temp_range): # from HW10. might be useful
        if temp_range == 'high':
            query = '''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7 FROM HIGH WHERE SPECIES_NAME="{}"'''.format(species_name)
        elif temp_range == 'low':
            query = '''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7 FROM LOW WHERE SPECIES_NAME="{}"'''.format(species_name)
        coeffs_list = cursor.execute(query).fetchall()
        coeffs = np.asarray(coeffs_list)
        return coeffs