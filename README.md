# Arnold-tongue-project
This repository includes codes used to simulated perception and perceptual learning of figure-ground segregation. 
generate_stimulus.m is written to build a 20*20 matrix of input contrasts.
generate_weights.m is written to build a weighting matrix for input contrasts, whigh takes into account the receptive fiold size of each grid vortex with respect to its eccentricity. 
VF2Cort.m converts visual field coordinates to cortical coordinates.
schematic_perception.m simulates the first session of experiment.
schematicNEW.m simulates all eight training sessions with network structural changes according to a hebbian learning rule.

