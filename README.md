# networkEv
## Purpose
This module allows the user to easily run an evolutionary algorithm to create mass-action chemical reaction networks with oscillatory behavior.

## How to use
### Installing Julia
This module uses the julia programming language. Installation instructions can be found [here](https://julialang.org/downloads/).

### Quick Start
1. Either fork or clone this module. Clone via the command ```git clone https://github.com/ltatka/networkEv```
2. Navigate to the directory where you have downloaded the networkEv repo.

#### Command Line
1. To install the networkEv module and all dependencies: ```julia --project=. -e 'using Pkg; Pkg.instantiate()'
2. To run evolution with default settings: ```julia --project=. evscript.jl <br>
This command also takes optional command line arguments: <br>
* ```--nbatches```: the number of batches to run (100 by default)
* ```--ngenerations```: the number of generations per batch (800 by default)
* ```--populationsize```: the number of networks in a population
* ```--pathtosettings```: path to a json file storing additional custom settings
* ```--outputpath```: path to a directory where evolution output will be written
* ```--seed```: specify a seed for the random number generator

#### Bash Script
```./quickstart.sh``` to run evolution with default settings. This script can also be edited to include the command line arguements in the previous section.

#### Write your own scripts
1. Use the run_evolution in a .jl file as follows <br>
``` using networkEv <br>
networkEv.run_evolution```

2. Command line arguments can be specified as keyword arguments in the ```run_evolution``` function.
For example ```networkEv.run_evolution(ngenerations=500)```


## Customizing Settings
Almost all settings can be customized. 

### Keyword Arguments
Some settings can be set from the julia command line. These are: <br>
```ngenerations```: The number of generations <br>
```nbatches```: The number of times to run the evolutionary algorithm <br>
```populationsize```: The number of reaction networks in a population <br>
```seed```: Set the random seed <br>
```pathtosettings```: A path to a JSON file where additional custom settings are specified <br>
```outputpath```: The path to which evolved reaction networks will be written <br>

Any of these keyword arguments can be specified when running evolution from the julia command line. For example:
```run_evolution(ngenerations=300, nbatches=100)```

### Additional Custom Settings
#### Specify Additional Settings in a JSON File
Additional settings can be specified in a JSON file. The path must then be input to the run_evolution function on the julia command line.
```run_evolution(pathtosettings="/home/name/path/to/your/settings.json")```

An example of a JSON file specifying custom settings: <br>
```
{"specieslist": ["A", "B", "C"],
"initialconditions": [5.0, 3.0, 6.0],
"reactionprobabilities": [0.2, 0.3, 0.3, 0.2],
"ngenerations": 500,
"populationsize": 200,
}
```
Any settings that are not specified in the JSON file will be set to the default value.

#### Description of Settings and Default Values
<table>
  <tr>
    <th>Setting Name</th>
    <th>Default Value</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>ngenerations</td>
    <td>400</td>
    <td>The number of generations</td>
  </tr>
  <tr>
    <td>populationsize</td>
    <td>100</td>
    <td>The number of reaction networks in the population</td>
  </tr>
  <tr>
    <td>nreactions</td>
    <td> 5 </td>
    <td>The number of reactions in each network at the beginning of evolution. This number may be different in final networks.</td>
  </tr>
  <tr>
    <td>target_num_species</td>
    <td>10</td>
    <td>The threshold to split reaction networks into species will be adjusted up and down to attempt to maintain this number of species. </td>
  </tr>
 IN PROGRESS


