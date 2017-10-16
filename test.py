import chemkin

def test_xmlErrors():
	try:
		chemkin.DataParser().parse_file('data/test_parse.xml') # no phase labels
		assert(False)
	except Exception as err:
		assert(type(err)==Exception)

def test_negA():
	try:
		chemikin.DataParser().parse_file("data/test_negA.xml")
		assert(False)
	except Exception as err:
		assert(type(err)==Exception)

def test_btype():
	try:
		chemikin.DataParser().parse_file("data/test_btype.xml")
		assert(False)
	except Exception as err:
		assert(type(err)==Exception)

def test_Atype():
	try:
		chemikin.DataParser().parse_file("data/test_Atype.xml")
		assert(False)
	except Exception as err:
		assert(type(err)==Exception)

def test_negk():
	try:
		chemikin.DataParser().parse_file("data/test_negk.xml")
		assert(False)
	except Exception as err:
		assert(type(err)==Exception)

def test_errorChem():
	try:
		chemikin.DataParser().parse_file("data/test_errorChem.xml")
		assert(False)
	except Exception as err:
		assert(type(err)==Exception)

def test_oneMoreReaction():
	try:
		chemikin.DataParser().parse_file("data/test_oneMoreReaction.xml")
	except Exception as err:
		assert(False)

def test_idCollision():
	try:
		chemikin.DataParser().parse_file("data/test_idCollision.xml")
	except Exception as err:
		assert(False)

