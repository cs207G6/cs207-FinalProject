import chemkin

def test_xmlErrors():
	try:
		chemkin.DataParser().parse_file('rxns.xml')
		chemkin.DataParser().parse_file('test_parse.xml')
	except Exception as err:
		assert(type(err)==Exception)

def test_zerocoeff():
	try:
		chemkin.DataParser().parse_file('rxns.xml')