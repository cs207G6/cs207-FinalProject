import os
import sqlite3

from . import nasa as n


class NASACoeffs:
    def __init__(self, database_file=None):
        if database_file is None:
            path = os.path.dirname(n.__file__)
            self.database_file = os.path.join(path, "database", "default.sqlite")
        else:
            self.database_file = database_file

    def create_db(self, filename):
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
        db = sqlite3.connect(self.database_file)
        # create a cursor object
        cursor = db.cursor()
        cursor.execute("DROP TABLE IF EXISTS high")
        cursor.execute("DROP TABLE IF EXISTS low")
        cursor.execute("DROP TABLE IF EXISTS all_temps")
        cursor.execute('''CREATE TABLE high (
                       SPECIES_NAME TEXT PRIMARY KEY NOT NULL, 
                       TLOW REAL, 
                       THIGH REAL, 
                       `COEFF_1` REAL,
                       `COEFF_2` REAL,
                       `COEFF_3` REAL,
                       `COEFF_4` REAL,
                       `COEFF_5` REAL,
                       `COEFF_6` REAL,
                       `COEFF_7` REAL)''')
        cursor.execute('''CREATE TABLE low (
                       SPECIES_NAME TEXT PRIMARY KEY NOT NULL, 
                       TLOW REAL, 
                       THIGH REAL, 
                       `COEFF_1` REAL,
                       `COEFF_2` REAL,
                       `COEFF_3` REAL,
                       `COEFF_4` REAL,
                       `COEFF_5` REAL,
                       `COEFF_6` REAL,
                       `COEFF_7` REAL)''')
        import xml.etree.ElementTree as et
        with open("./test.xml") as f:
            root = et.fromstring(f.read())
        speciesData = root.find("speciesData")
        high_data = []
        low_data = []
        for species in speciesData.findall("species"):
            sname = species.attrib['name']
            NASAs = species.find('thermal').findall('NASA')
            assert (len(NASAs) == 2)
            max0 = NASAs[0].attrib['Tmax']
            max1 = NASAs[1].attrib['Tmax']
            if max0 > max1:
                high = NASAs[0]
                low = NASAs[1]
            else:
                high = NASAs[1]
                low = NASAs[0]
            get_data = lambda d: tuple(
                [sname, d.attrib['Tmin'], d.attrib['Tmax']] + [float(e) for e in d.find('floatArray').text.split(',')])
            high_data.append(get_data(high))
            low_data.append(get_data(low))

        cursor.executemany('''INSERT INTO high 
                          (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7)
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', high_data)
        cursor.executemany('''INSERT INTO low 
                          (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7)
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', low_data)
        db.commit()
        cursor.execute('''SELECT * FROM high''').fetchall()
        db.close()
        cursor.close()

    def get_coeffs(self, species_name, temp_range):
        db = sqlite3.connect(self.database_file)
        cursor = db.cursor()
        # get list of species
        query = '''SELECT SPECIES_NAME FROM LOW'''
        species_list = cursor.execute(query).fetchall()
        species_list = [s[0].strip(',') for s in species_list]
        # check that specie_name inputed it correct type of species
        if species_name not in species_list:
            raise ValueError("Input must be 'O','O2','H','H2','OH','H2O','HO2','H2O2' for species_name")
        # get the coeffs based on temp range
        if temp_range == 'low':
            query = '''SELECT COEFF_1,COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7,TLOW,THIGH FROM LOW WHERE SPECIES_NAME = ? '''
        elif temp_range == 'high':
            query = '''SELECT COEFF_1,COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7,TLOW,THIGH FROM HIGH WHERE SPECIES_NAME = ? '''
        else:  # make sure temp range is correctly inputed
            raise ValueError("Must input 'low' or 'high' for temp_range")
        results = cursor.execute(query, (species_name,)).fetchall()
        if len(results) != 1:
            raise KeyError("NASA coefficients for {} is not defined".format(species_name))
        result = [float(i) for i in results[0]]
        coeffs = result[:7]
        tmin = result[7]
        tmax = result[8]
        cursor.close()
        db.close()
        return coeffs, tmin, tmax
