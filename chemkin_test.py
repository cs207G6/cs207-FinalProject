import chemkin

def test_xmlErrors():
	try:
		chemkin.DataParser().parse_file('data/test_parse.xml') # no phase labels
		assert(False)
	except Exception as err:
		assert(True)

def test_negA():
	try:
		chemkin.DataParser().parse_file("data/test_negA.xml")
		assert(False)
	except Exception as err:
		assert(True)

def test_btype():
	try:
		chemkin.DataParser().parse_file("data/test_btype.xml")
		assert(False)
	except Exception as err:
		assert(True)

def test_Atype():
	try:
		chemkin.DataParser().parse_file("data/test_Atype.xml")
		assert(False)
	except Exception as err:
		assert(True)

def test_negk():
	try:
		chemkin.DataParser().parse_file("data/test_negk.xml")
		assert(False)
	except Exception as err:
		assert(True)

def test_errorChem():
	try:
		chemkin.DataParser().parse_file("data/test_errorChem.xml")
		assert(False)
	except Exception as err:
		assert(True)

def test_oneMoreReaction():
	try:
		chemkin.DataParser().parse_file("data/test_oneMoreReaction.xml")
	except Exception as err:
		assert(False)

def test_length():
    assert(len(chemkin.DataParser().parse_file('data/rxns.xml'))==3)

def test_idCollision():
	try:
		chemkin.DataParser().parse_file("data/test_idCollision.xml")
		assert(False)
	except Exception as err:
		assert(True)
        
