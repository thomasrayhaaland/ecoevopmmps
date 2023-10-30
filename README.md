## README for eco_evo_pmmps.R ##
Code accompanying Haaland et al.: Eco-evolutinary dynamics of partially migratory metapopulations

# Simulation program #
The main simulation program begins on line 1087. It requires the preamble (chunk line 1-63) to be run before executing.
It includes sensible defaults and so should run smoothly immediately.
A list of arguments and outputs are provided above (from line 1028).

# Simulation data #
Simulation output must be saved either locally or to disk - see example on line 1460 & 1461.
Before being used for plotting, simulation output needs to be converted to a stats file.
This is done with the popstats function. Run the chunk "Functions for wrangling and plotting results" (line 65-870) first.
Then take simulation output from local environment or disk (line 1482) and save the statsfile locally or to disk (line 1483).

# Plotting #
Plotting occurs using popplots function - see line 171 and arguments list above.
Line 1492-1826 gives examples of code for plotting publication figures.

# Individual-level databases #
Deatiled yearly data is collated into data file using functions create.alldb (line 880) and polish.alldb (line 918). 
This must be saved as a separate alldb file before using for analyses - see example line 1838.
Evolutionary dynamics data file are created from alldb using function create.evoldyn - line 1855.
See further examples of plotting after this.

