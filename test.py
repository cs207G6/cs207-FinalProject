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

def test_negk():
	try:
		chemikin.DataParser().parse_file("data/test_negk.xml")
		assert(False)
	except Exception as err:
		assert(type(err)==Exception)

