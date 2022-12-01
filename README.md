# SqueezePIR

## Environment
SqueezePIR is built on SEAL 3.7.0, cmake version is 3.10


## Run the code
```
cmake . -DCMAKE_PREFIX_PATH=~/direct_to_Microsoft_SEAL_build
make
./bin/SqueezePIR
```

## Parameter description
In this work, the degree of poly in CKKS is 8192
NUM_OBJ number of rows for the database (2^n)
NUM_COL number of columns for the database      
RANK    rank of SqueezePIR          
IDX     entry ID that the client wants to retrieve
