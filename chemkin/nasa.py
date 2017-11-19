import sqlite3

class NASACoeffs:
    def __init__(self, database_file = None):
        if database_file is None:
            self.database_file = "default_database"
        else:
            self.database_file = database_file
    
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
        db = sqlite3.connect(self.database_file)
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
        cursor.executemany('''INSERT INTO LOW
                    (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', lows_to_insert)
        cursor.executemany('''INSERT INTO HIGH
                    (SPECIES_NAME, TLOW, THIGH, COEFF_1, COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', high_to_insert)
    
        db.commit()
        db.close()

    def get_coeffs(species_name, temp_range):
        #get list of species
        query = '''SELECT SPECIES_NAME FROM LOW'''
        species_list = cursor.execute(query).fetchall()
        species_list = [s[0].strip(',') for s in species_list]
        #check that specie_name inputed it correct type of species
        if species_name not in species_list:
            raise ValueError("Input must be 'O','O2','H','H2','OH','H2O','HO2','H2O2' for species_name")
        #get the coeffs based on temp range
        if temp_range == 'low':
            query = '''SELECT COEFF_1,COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7 FROM LOW WHERE SPECIES_NAME = ? '''
            coeffs = cursor.execute(query,(species_name,)).fetchall()
        elif temp_range == 'high': 
            query = '''SELECT COEFF_1,COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7 FROM HIGH WHERE SPECIES_NAME = ? '''
            coeffs = cursor.execute(query,(species_name,)).fetchall()
        else: #make sure temp range is correctly inputed
            raise ValueError("Must input 'low' or 'high' for temp_range")
        return coeffs
    
    def get_tmin(species_name, temp_range):
         #get the tmin based on temp range
        if temp_range == 'low':
            query = '''SELECT TLOW FROM LOW WHERE SPECIES_NAME = ? '''
            tmin = cursor.execute(query,(species_name,)).fetchall()
        elif temp_range == 'high': 
            query = '''SELECT TLOW FROM HIGH WHERE SPECIES_NAME = ? '''
            tmin = cursor.execute(query,(species_name,)).fetchall()
        else: #make sure temp range is correctly inputed
            raise ValueError("Must input 'low' or 'high' for temp_range")
        return tmin
    
    def get_tmax(species_name, temp_range):
         #get the tmax based on temp range
        if temp_range == 'low':
            query = '''SELECT THIGH FROM LOW WHERE SPECIES_NAME = ? '''
            tmax = cursor.execute(query,(species_name,)).fetchall()
        elif temp_range == 'high': 
            query = '''SELECT THIGH FROM HIGH WHERE SPECIES_NAME = ? '''
            tmax = cursor.execute(query,(species_name,)).fetchall()
        else: #make sure temp range is correctly inputed
            raise ValueError("Must input 'low' or 'high' for temp_range")
        return tmax
