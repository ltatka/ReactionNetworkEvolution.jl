#!/bin/bash
counter=1 

while [ $counter -le 7 ] 

do  

    julia --project=. -e 'using Pkg; Pkg.instantiate()'  

    julia --project=. run_evolution.jl --pathtosettings="/home/hellsbells/Desktop/evolution_output/no_speciation_new/settings.json" --outputpath="/home/hellsbells/Desktop/evolution_output/no_speciation_new/"

    ((counter++)) 

done 
