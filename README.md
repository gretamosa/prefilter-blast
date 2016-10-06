Prefilter BLAST
================
[Email](mailto:german.retamosa@uam.es)

Prefilter BLAST is a novel parallel prefiltering model written for NVidia GPU. 
Although it could be used with other algorithms, it has been specially developed and optimized for NCBI BLAST.

![](https://images.duckduckgo.com/iu/?u=http%3A%2F%2Fstatic.guim.co.uk%2Fsys-images%2FGuardian%2FPix%2Fcartoons%2F2013%2F9%2F9%2F1378724614866%2Fgattaca-genetic-sequencin-008.jpg&f=1)

## Features
### Protein Filtering
- Stable version

### Nucleotide Filtering
- Working in progress (alpha development version)

### Metagenome
- Working in progress (alpha development version)

## Requirements
Prefilter BLAST has been developed for NVidia GPU and its programming language CUDA.
For that reason, CUDA and its dependecies must be installed. For further information, go to [NVidia Developer Portal](http://docs.nvidia.com/cuda/cuda-quick-start-guide/index.html#linux-x86_64)

## Dependencies
- CUDA Runtime & Developer Toolkit (tested with CUDA 2.0) 
- GCC compiler and development dependencies (glibc-dev)

## Parameters
- -p: NCBI BLAST type (blastp for proteins, blastn for nucleotides and blastx for metagenomes).
- -d: Database file (i.e. non-redundant NR) 
- -i: Query file (i.e. AMA protein family)
- -b: Filter mode (1 for different alignments, 2 for similar alignments; default 1)
- -o: Output FASTA file
- -f: Prefilter BLAST similarity threshold (default 0.5)
- -t: NCBI BLAST threshold (default 11)
- -c: GPU card number (specially for multi-GPU cards; default 0)
- -m: NCBI BLAST substitution matrix (i.e. BLOSUM62)
- -w: NCBI BLAST word size (3 for proteins and 11 for nucleotides)

## Compilation
```bash
make
```

## Execution
```bash
naudit_filter -p blastp -d db/nr.fasta.fdb -i db/Amaaa.fasta -b 2 -o salida -a 6 -t 0.7
```

## Contribute
If you have any idea for an improvement or found a bug do not hesitate to open an issue.

## License
Prefilter BLAST is distributed under MIT License.
