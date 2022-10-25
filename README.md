# Code for paper: Towards the Scalability of Distributed Network Emulations
## Build
```
mkdir build
cd build
cmake .. && make 
```

## Run
```
./tbs <topology file> --k=<number of parts> --preconfiguration=esocial
```

## Convert to MaxiNet nodemapping
Run
```
python convert.py
```
to generate a python dictionary object and pass as the argument of `nodemapping` on the creation of  `Experiment` object in MaxiNet.

## Disclaimer
Some codes here are developed upon [KaHIP](https://github.com/KaHIP/KaHIP)
