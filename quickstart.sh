#!/bin/bash
counter=1 

while [ $counter -le 5 ] 

do  

    julia --project=. -e 'using Pkg; Pkg.instantiate()'  

    julia --project=. evscript.jl

    ((counter++)) 

done 
