UI, visualization model, expand NASA coefficients


## UI
- Upload file (with one or multiple reactions)
	- .xml
	- .txt
	- .HDF5
	
- Manual Inputs
	- Temperature
	- Species
	- Concentrations
	- Whether reversible
	- Elementary reactions
		- Reaction Coefficient
			- A,b,E OR
			- Direct k
		- Reactants
			- species name (choose from inputed species)
			-  stoichiometric coefficients
		- Products
			- species name (choose from inputed species)
			-  stoichiometric coefficients

- Plot choices
	- whether to plot or not 
	- Given temperature range (optional)
	- only plot certain species 


## Python module
### Main
- call data processing function
	- Arguments: user inputs
	- Return:
		- .xml
		- Indicator variables about plots
 
- call visualization function
	- Arguments: output from data processing
	- Return: plots
 
### Data processing
- Convert the user inputs from webUI to .xml (with normal format)
- Return: .xml file and indicator variables (e.g. is there a given range of temp; is there certain wanted species)


### Visualization
- uploaded a file
	- read the temperature in the file
- Given temperature range
- Otherwise
	-  pick the nasa high/low range

## Expand NASA Coeffs
- DONE!!!?
