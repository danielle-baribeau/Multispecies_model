### MANAGEMENT SCENARIO FUNCTION ###

#This function is called in trophic_level_RUNME_mgmt, and fills in one of a set of default management scenarios

# DB Notes
#1: type:     #integer representing which management scenario type you want to run
  #(how u will change with biomass between biomass thresholds)
  #options:
    #1 - u is constant regardless of biomass
    #2 - u at constant minimum below low threshold, increases linearly between thresholds, 
      #then at constant max above high threshold (DFO standard protocol)
    #3 - u at constant minimum below low threshold, increases linearly beyond low threshold

mgmt.scen.setup <- function()
{}