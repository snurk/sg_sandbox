# String Graph Sandbox
Semi-automated prototype pipeline for assembling (pseudo-)haploid genomes, sex chromosomes or highly inbred genomes from HiFi (and ONT data).

This repo brings together experimental procedures for:
* construction of string graph based on HiFi reads
* pruning of the graph
* manual-assisted ONT-based repeat resolution
* generation of consensus sequences based on the specified graph paths

An early version of this pipeline has been used in assembling a telomere-to-telomere draft of CHM13 human genome.

## Requirements
Requirements might differ depending on which part(s) of the pipeline you are planning to use, but you might probably need:

* Modern C++ compiler
* Seqtk, minimap2
* `Python3` environment with `biopython`, `parasail-python` and `pysam`
* environment for compiling GraphAligner see [here](https://github.com/maickrau/GraphAligner)

## Compilation
Call `build.sh` within the environment that can compile GraphAligner.

## Building and pruning the graph
