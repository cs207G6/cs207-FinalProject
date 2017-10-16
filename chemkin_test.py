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

def test_progRateWrongDimension():
	cd = chemkin.DataParser().parse_file('data/rxns.xml')
	try:
		cd.get_progress_rate([1,2,3,4,5],100)
		assert(False)
	except Exception as err:
		assert(True)
        
def test_progRateWrongSp():
	try:
		chemkin.DataParser().parse_file('data/test_wrongSp.xml')
		assert(False)
	except Exception as err:
		assert(True)
        

def test_progRateWrongSp2():
	try:
		chemkin.DataParser().parse_file('data/test_wrongSp2.xml')
		assert(False)
	except Exception as err:
		assert(True)
        

def test_progRateNonIrr():
	cd = chemkin.DataParser().parse_file('data/test_nonirr.xml')
	try:
		cd.get_progress_rate([1,2,3,4,5,6],100)
		assert(False)
	except Exception as err:
		assert(True)
        
def test_progRateNonEle():
	cd = chemkin.DataParser().parse_file('data/test_nonele.xml')
	try:
		cd.get_progress_rate([1,2,3,4,5,6],100)
		assert(False)
	except Exception as err:
		assert(True)
        