
# coding: utf-8

# In[67]:


class ReactionData():
    def __init__(self, id, species, reactions):
        self.id = id
        self.reactions = reactions
        self.species = species

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
    
    def get_K(self):
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
        if k < 0:
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
    
    def parse_rate_coeff(self, root):
        root = root.getchildren()[0]
        if root.tag == "Arrhenius":
            A = float(root.find('A').text)
            E = float(root.find('E').text)
            return Arrhenius(a=A, E=E)
        elif root.tag == "modifiedArrhenius":
            A = float(root.find('A').text)
            b = float(root.find('b').text)
            E = float(root.find('E').text)
            return ModifiedArrhenius(a=A, b=b, E=E)
        elif root.tag == "Constant":
            k = float(root.find('k').text)
            return Constant(k)
        else:
            raise ValueError("[TODO]")
    
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
    
    def parse_reaction(self, root):
        rate_coeff = self.parse_rate_coeff(root.find('rateCoeff'))
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
            reaction = self.parse_reaction(r)
            reactions.append(reaction)
        return ReactionData(id, species, reactions)


# In[69]:


reactionData = DataParser().parse_file("data/rxns.xml")


# In[72]:


reaction0 = reactionData.reactions[0]

