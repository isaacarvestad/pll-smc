# pll-smc
Sequential Monte Carlo phylogenetic tree inference using
the [PosetSMC](https://www.stat.ubc.ca/~bouchard/PosetSMC/) algorithm with
likelihood computations performed by
the [Phylogenetic Likelihood Library](https://github.com/xflouris/libpll).

## Getting started
Instructions assume a Unix system with Bash or a similar shell.

### Prerequisites
Before building, make sure the following is installed on your system:
- [libpll](https://github.com/xflouris/libpll)
- [CMake](https://cmake.org)

### Building
When in the root project directory:

``` bash
mkdir build
cd build
cmake ..
make
```

This will build a runnable binary at `[...]/pll-smc/build/app/pll-smc`.

### Usage
To run pll-smc, provide a path to a set of sequences in the Fasta format.

``` bash
# Assuming inside 'build' directory
./app/pll-smc path/to/sequences.fasta
```

Optionally, the number of particles can be specified. If a value isn't provided
default is `1000`.

``` bash
# Assuming inside 'build' directory
./app/pll-smc path/to/sequences.fasta 500
```

Once the tree distribution has been inferred it will be written to
stdout. Progress information is written to stderr continuously during
execution. To save the tree distribution we can redirect it to a file.

``` bash
# Assuming inside 'build' directory
./app/pll-smc path/to/sequences.fasta > output.txt
```

## Author
- Isaac Arvestad

## Contributors
- Henrik Lagebrand

## License
See [license](LICENSE).
