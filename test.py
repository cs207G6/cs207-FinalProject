import chemkin

def test_xmlErrors():
	try:
		#chemkin.DataParser().parse_file('data/rxns.xml')
		chemkin.DataParser().parse_file('data/test_parse.xml') # no phase labels
	except Exception as err:
		assert(type(err)==Exception)

def test_negA():
	try:
		#reactionData0 = chemikin.DataParser().parse_file("data/rxns.xml")
		chemikin.DataParser().parse_file("data/test_negA.xml")
	except Exception as err:
		assert(type(err)==Exception)

def test_btype():
	try:
		#reactionData0 = chemikin.DataParser().parse_file("data/rxns.xml")
		chemikin.DataParser().parse_file("data/test_btype.xml")
	except Exception as err:
		assert(type(err)==Exception)

