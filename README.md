## Table of Contents

- [Users' Guide](#uguide)
  - [Dependecies](#install)
  - [General usage](#general)
  - [Map short accurate genomic reads](#short-genomic)
- [Algorithm overview](#algo)
- [Limitations](#limit)

## <a name="uguide"></a>Users' Guide

MapTie is a versatile sequence alignment program that aligns DNA sequences against a large reference genome. 


### <a name="install"></a>Dependecies 

```sh
dependencies:
  - python=3.7
  - scipy
  - numpy
  - matplotlib
  - numba
```
### <a name="general"></a>Configuration

MapTie requires a configuration file conformant with this standard:

```sh
{
    "input":{
        "reference_abs_path": "/Users/sat/Desktop/ETH_PROJECTS/comp_bio/team02/data_small/genome.chr22.5K.fa", # Reference genome absolute path
        "reference_name": "chr22_tiny", 
        "read1": "/Users/sat/Desktop/ETH_PROJECTS/comp_bio/team02/data_small/output_tiny_30xCov1.fq", # First read file absolute path
        "read2": "/Users/sat/Desktop/ETH_PROJECTS/comp_bio/team02/data_small/output_tiny_30xCov2.fq"  # Second read file absolute path
    },
    "output":{
        "cache_dir": "cache", 
        "output_dir": "results"
    },
    "execution":{
        "n_jobs":2  # Maximum number of concurrently running workers
    }
}

```

### <a name="short-genomic"></a>Map short accurate genomic reads

Once the config.json file is configured, mapping the reads is as simple as:

```sh
python3 maptie.py     # paired-end alignment
```
When two read files are specified, MapTie reads from each file in turn and
merge them into an interleaved stream internally. Two reads are considered to
be paired if they are adjacent in the input stream and have the same id. 

MapTie does not work with unpaired reads.




## <a name="algo"></a>Algorithm overview

Maybe some details about the algorithm here.




## <a name="limit"></a>Limitations

Maybe some details about "weak-spots" of the mapper.



