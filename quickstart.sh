#!/bin/bash
counter=1 

while [ $counter -le 5 ] 

do  

    julia --project=. -e 'using Pkg; Pkg.instantiate()'  

    julia --project=. evscript.jl --ngenerations 1 --nbatches 1

    ((counter++)) 

done 
