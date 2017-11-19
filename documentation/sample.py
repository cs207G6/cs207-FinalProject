import chemkin

data_parser = chemkin.DataParser() # create a data parser class
reaction_data = data_parser.parse_file("data/rxns.xml") # parse the data file and return an instance of ReactionData class
progress_rates = reaction_data.get_progress_rate([1,2,3,4,5,6],100)
print(progress_rates)
reaction_rates = reaction_data.get_reaction_rate(progress_rates)
print(reaction_rates)
ks = reaction_data.get_k(20)
print(ks)