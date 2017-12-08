import chemkin.parser
import chemkin.nasa

nasa = chemkin.nasa.NASACoeffs()
data_parser = chemkin.parser.DataParser()

reaction_data = data_parser.parse_file("chemkin/example_data/example_thermo.xml", nasa)

# Get Progress Rate
progress_rates = reaction_data.get_progress_rate([1, 2, 3, 4, 5, 6], 100)
print(progress_rates)

# Get Reaction Rate
reaction_rates = reaction_data.get_reaction_rate(progress_rates)
print(reaction_rates)

# Get Reaction Rate Coefficients
ks = reaction_data.get_k(20)
print(ks)
