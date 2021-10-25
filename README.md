# String Graph Sandbox
Semi-automated prototype pipeline for assembling (pseudo-)haploid genomes, sex chromosomes or highly inbred genomes from HiFi (and ONT data).

In all honesty, the pipeline is so experimental that this documentation might not be enough for effectively assembling anything without contacting at least one of the Sergeys :)

This repo brings together experimental procedures for:
* construction of string graph based on HiFi reads
* pruning of the graph
* manual-assisted ONT-based repeat resolution
* generation of consensus sequences based on the specified graph paths

An early version of this pipeline has been used in assembling a telomere-to-telomere draft of CHM13 human genome.

## Requirements
Requirements differ depending on which part(s) of the pipeline you are planning to use, but you might need:

* Modern C++ compiler
* `GraphAligner` (v1.0.13+), `seqtk`, `minimap2`, `samtools` and `bedtools` in PATH
* `Python3` environment with `snakemake`, `biopython`, `parasail-python` and `pysam` (see also `requirements.txt`)
* environment for compiling GraphAligner see [here](https://github.com/maickrau/GraphAligner)

## Compilation
Call `build.sh`.

## Building the graph

1. First, you might consider selecting a maximum read-length to trim longer reads to which mitigates inherent issues of string graphs built from reads of varying lengths.
1. Then, tweak `config.yaml` from `src/pipe/`. You will need to edit it to set the S_MAX_RL to the trim length you used above adjusted for homopolymer compression (something like `trim length` / 1.4 * 1.1). You may also need to adjust the S_EXPECTED_COV to (probably slightly below) the (haploid if you have a diploid non-inbred sample) coverage that you expect to have for 'unique' regions of your chromosome(s).
1. `canu.sh` from `src/canu_launch` might also need tweaking. First, for performance reasons if you are working with something much larger than the human genome. Second, to adjust grid options to match your cluster requirements.
1. Depending on your cluster configuration you might also need to modify submission commands in `src/canu_launch/on_init_complete.sh` and `src/canu_launch/on_primary_complete.sh`, as well as `src/pipe/cluster.json` and `src/pipe/launch.sh` for Snakemake to correctly submit the graph construction and processing jobs.
1. Run `src/canu_launch/master_trim.sh` script, which will 
    1. Run read trimming
    2. Run two iterations of read/overlap processing from HiCanu
    3. Generate `updated_raw.seqStore` necessary for final consensus generation
    4. Meyers' string graph construction using modified `miniasm` code
    5. Pruning the graph using custom procedures from the `gfacpp` repo (see `src/pipe/simplif.sh`)

When done, you will have a `simplified.gfa` (as well as `simplified.noseq.gfa` and `simplified.nodes.fasta` in your specified output folder, which is the primary outputs of this part of the pipeline.
Note that all the sequences are obtained from homopolymer-compressed version of the reads. 
So unitigs themselves will be (mostly) homopolymer-free (neither their sequences nor lengths directly translate to those of underlying genomic sequences).
Another output of notice is `resolved_mapping.txt`, providing information about the non-contained-read-level layout of individual graph unitigs.

## Visualizing the graph
To visualize the graph you can use the amazing [Bandage](https://github.com/rrwick/Bandage) tool .
If you are planning to use ONT-based resolution you will need to curate a list of `unique` nodes that have genomic multiplicity 1.
Node lengths and coverage estimates (included in `simplified.gfa`) can help with getting an approximation, but you will probably need to vet the list after exploring the graph in Bandage.

## 'Semi-automated' ONT-based resolution
This repo includes the original ONT-based repeat resolution strategy developed by Mikko during our CHM13 finishing workshop, slightly modified to be used with string graphs.
To use it you will need to first 'compress' the homopolymer runs in your ONT reads (you can use the `./src/contig_processing/homopolymer_compress.py`, but it is a bit slow) and probably length-filter them to speed things up.

Briefly it will:
1. Align homopolymer-compressed ONT reads to the graph using GraphAligner
2. Process alignment paths to disambiguate and trim unreliable path ends
3. For each unique node identify unique nodes following/preceding it in the alignment paths
4. Filter out 'bridges' that have a better supported alternative (at least twice as many reads supporting the connection)
5. Look at connected components of the graph built on the sides of unique nodes with edges corresponding to reachability via non-unique nodes
6. If all the unique node sides have at least one 'bridge' -- resolve component. Namely,
7. Connect unique node sides according to alignment paths forming the 'bridge' (making copies of non-unique nodes as necessary)
8. If there are multiple different paths connecting same unique nodes -- pick the one which is 2 times better supported than any alternative. If such path doesn't exist -- introduce multiple connecting paths.

## 'Manual' resolution and graph tweaking
At the moment we don't provide the documentation for that part.
And hope that nobody will need to use it now that a new generation of automated solutions is being developed.
If you really want to go along a similar path to one we've been walking -- contact us.

## Re-compacting the graph and forming layouts
If you have done ONT resolution and/or manual editing of the GFA in Bandage you might want to 're-compact' the graph, e.g. to facilitate further layout picking or another round of ONT resolution.
This can be done with the `src/process_simplified.sh`.

Now you will need to form the `layout.txt` file(s) specifying the paths for which you want to compute the consensus sequences.

Each line of the `layout.txt` file should have the following format: `<name> <node1>[+,-],<node2>[+,-],...`

**NB:** Bangage's *Output/Specify exact path* functionality is very convenient for selecting multi-node paths.

## Consensus and post-processing
### Consensus
Consensus for the identified layouts can be obtained by running `src/consensus/launch.sh`, providing the working folder with `layout.txt` file as one of the parameters. The resulting output will be in `cns.renamed.fasta` file in the that folder.

Provided `src/consensus/config.yaml` should not need any modification if you were building the graph via the `src/canu_launch/master_trim.sh`.

You might need to edit `src/consensus/launch.sh` and/or `src/consensus/cluster.json` to adjust it to your cluster configuration and submission system.

### Layout splitting
Sometimes, due to the limitations of the consensus module, the layouts need to be additionally split on repetitive nodes formed by few non-contained reads. If the 'layouReads' binary failed then you most likely hit this issue.
Layout splitting can be done mostly automatically using the `consensus/backbone_analysis.py` script.

### Patching
Sometimes you might want to use a different 'reference' assembly to patch the gaps between your contigs.
In case of our work, we have been using Flye assemblies of the ONT reads to patch the gaps originating from the 'GA'-microsatellite dropout issue of HiFi sequencing.
This can be done with the `src/contig_processing/master_patch.sh` script, which will align your contigs to the 'reference' assembly with minimap2 and will try to use the alignments to then fill in the missing regions with the reference regions. 

**NB.** After running the script check the `patch.bed` file in the output folder to see if any gaps are unfilled.
Basically, your contig records should alternate with reference contig records. For example:
```
chrX.01    0    191138    .    0    +
contig_4387    679516    709990    .    0    -
chrX.02    0    57646    .    0    +
```
is OK, BUT

```
chrX.01    0    191138    .    0    +
chrX.02    0    57646    .    0    +
```
is NOT!

Unfortunately, in the latter case (which can happen if the contigs actually overlap; or the gap was too big; or alignments could not be identified) no error message will be generated and the contigs will be erroneously concatenated.

### Stitching
If during the patching step (or if you had to 'split' your layout) you realized that your resulting contigs should actually overlap you can 'stitch' them using `sg_sandbox/src/consensus/join_ctgs.py` script, which will try to identify overlap between the contigs and then merge them.
Among other parameters you will need to provide it with a stitch 'plan' to indicate the order of the contigs you would like to try concatenating.
The format is the same as for the layout.txt files, e.g. one or more lines of the format: `chrX.1 chrX.1.1+,chrX.1.2+,...`.
Confirm that the identified overlaps correspond to what you expected.

## Example
You can try running 24X coverage _E.coli_ dataset available [here](https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz).
```
src/canu_launch/master_trim.sh ~/git/sg_sandbox/src/pipe/config.yaml pipeline_test 23000 subset24x.fastq.gz
#Check that pipeline_test/simplified.noseq.gfa exists, non-empty and contains a single S-line. Copy segment_name.
mkdir -p pipeline_test/consensus/chr
cd pipeline_test/consensus/chr
echo 'chr <segment_name>+' > layout.txt
src/consensus/launch.sh CHR ~/git/sg_sandbox/src/consensus/config.yaml .
```
