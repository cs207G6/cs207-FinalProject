import chemkin

reaction_data = chemkin.DataParser().parse_file("data/rxns.xml")
progress_rates = reaction_data.get_progress_rate([1,2,3,4,5,6],100)
print(progress_rates)
reaction_rates = reaction_data.get_reaction_rate(progress_rates)
print(reaction_rates)
