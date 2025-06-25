! Warning This is an on-going project with incomplete analyses

# MANERRFC

Processing of Mission Aransas NERR and SWMP data for analysis investigating microplankton trophic roles across regime cycles.

# Organization

Scripts follow a category-numeric labeling scheme. Most data processing is all through R. 'pipe' scripts imports external data and formats them for local use. Scripts which conduct analyses are not prefixed, but follow the same numbering convention. Generated data files are stored as RDS, with the numeric code referencing the script which generated those data. External data which are available locally are denoted as 'xx-'. If a script requires external data which are not available, it should indicate in the top where to access those data. For direct questions or access, contact alex barth. 

Scripts which generate figures are denoted as 'fig' and will place items into the output folder. Stan code used for model building is located in a separate directory from R code.
