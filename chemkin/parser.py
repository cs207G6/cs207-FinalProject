from chemkin.reaction import *
from chemkin.rate_coeff import *


class DataParser:
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
                R = root.find('R')
                return Arrhenius(a=A, E=E, R=float(R.text))
            return Arrhenius(a=A, E=E)
        elif root.tag == "modifiedArrhenius":
            A = float(root.find('A').text)
            b = float(root.find('b').text)
            E = float(root.find('E').text)
            if root.find('R') is not None:
                R = root.find('R')
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
            result[parts[0]] = float(parts[1])
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
        return Reaction(id, reversible=reversible, type_=type_, reactants=reactants, products=products,
                        rate_coeff=rate_coeff)

    def parse_file(self, filename, nasa):
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
        for (i, r) in enumerate(tree.find('reactionData').findall('reaction')):
            reaction = self._parse_reaction(r)
            reactions.append(reaction)
        return ReactionData(id, species, reactions, nasa)
