import chemkin

def test_xmlErrors():
	try:
		chemkin.DataParser().parse_file('data/rxns.xml')
		chemkin.DataParser().parse_file('data/test_parse.xml')
	except Exception as err:
		assert(type(err)==Exception)

def test_zerocoeff():
	try:
		reactionData = chemikin.DataParser().parse_file("data/rxns.xml")
		reactionData
