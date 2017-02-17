#!/usr/bin/env python2.7
"""
vg_config.py: Default configuration values all here (and only here), as well as logic
for reading and generating config files.

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import getpass
import pdb
import textwrap
import yaml
from toil_vg.vg_common import require

default_config = textwrap.dedent("""
# Toil VG Pipeline configuration file (created by toil-vg generate-config)
# This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
# Edit the values in the configuration file and then rerun the pipeline: "toil-vg run"
# 
# URLs can take the form: "/", "s3://"
# Local inputs follow the URL convention: "/full/path/to/input.txt"
# S3 URLs follow the convention: "s3://bucket/directory/file.txt"
#
# Comments (beginning with #) do not need to be removed. 
# Command-line options take priority over parameters in this file.  
######################################################################################################################

###########################################
### Toil resource tuning                ###

# These parameters must be adjusted based on data and cluster size
# when running on anything other than single-machine mode

# TODO: Reduce number of parameters here.  Seems fine grained, especially for disk/mem
# option to spin off config files for small/medium/large datasets?   

# The following parameters assign resources to small helper jobs that typically don't do 
    # do any computing outside of toil overhead.  Generally do not need to be changed. 
misc-cores: 1
misc-mem: '1G'
misc-disk: '1G'

# Resources allotted for xg indexing.
# this stage generally cannot take advantage of more than one thread
# For whole genome, suggest 200G mem and disk
xg-index-cores: 1
xg-index-mem: '4G'
xg-index-disk: '2G'

# Resources allotted for gcsa pruning.  Note that the vg mod commands used in
# this stage generally cannot take advantage of more than one thread
# For whole genome, suggest 60G mem and disk
prune-cores: 1
prune-mem: '4G'
prune-disk: '2G'

# Resources allotted for gcsa kmers.  
# For whole genome, suggest 70G mem and disk
kmers-cores: 1
kmers-mem: '4G'
kmers-disk: '2G'

# Resources allotted gcsa indexing
# For whole genome, suggest 200G mem and 3T disk
gcsa-index-cores: 1
gcsa-index-mem: '4G'
gcsa-index-disk: '2G'

# Resources for fastq splitting and gam merging
# Important to assign as many cores as possible here for large fastq inputs
fq-split-cores: 1
fq-split-mem: '4G'
fq-split-disk: '2G'

# Number of threads to use for Rocksdb GAM indexing
# Generally, this should be kept low as speedup drops off radically 
# after a few threads.
gam-index-cores: 1

# Resources for *each* vg map job
# the number of vg map jobs is controlled by reads-per-chunk (below)
# For whole genome, suggest 100G memory and disk for whole genome
alignment-cores: 1
alignment-mem: '4G'
alignment-disk: '2G'

# Resources for chunking up a graph/gam for calling (and merging)
# typically take xg for whoe grpah, and gam for a chromosome,
# and split up into chunks of call-chunk-size (below)
# For whole genome, suggest 20G memory and 100G disk
call-chunk-cores: 1
call-chunk-mem: '4G'
call-chunk-disk: '2G'

# Resources for calling each chunk (currently includes pileup/call/genotype)
# For whole genome, suggest 20G memory and 20G disk
calling-cores: 1
calling-mem: '4G'
calling-disk: '2G'

# Resources for vcfeval
vcfeval-cores: 1
vcfeval-mem: '4G'
vcfeval-disk: '2G'

###########################################
### Arguments Shared Between Components ###
# Use output store instead of toil for all intermediate files (use only for debugging)
force-outstore: False

#############################
### Docker Tool Arguments ###
# Do not use docker for any commands
no-docker: False

## Docker Tool List ##
##   Each tool is specified as a list where the first element is the docker image URL,
##   and the second element indicates if the docker image has an entrypoint or not
##   If left blank or commented, then the tool will be run directly from the command line instead
##   of through docker. no-docker (above) overrides all these options. 

# Docker container to use for vg
vg-docker: ['quay.io/glennhickey/vg:v1.4.0-2213-gdce81f8', True]

# Docker container to use for bcftools
bcftools-docker: ['quay.io/cmarkello/bcftools', False]

# Docker container to use for tabix
tabix-docker: ['quay.io/cmarkello/htslib:latest', False]

# Docker container to use for jq
jq-docker: ['devorbitus/ubuntu-bash-jq-curl', False]

# Docker container to use for rtg
rtg-docker: ['realtimegenomics/rtg-tools:3.7.1', True]

# Docker container to use for pigz
pigz-docker: ['quay.io/glennhickey/pigz:latest', True]

##########################
### vg_index Arguments ###

# Name of index output files.  ex <name>.xg, <name>.gcsa etc. 
index-name: 'genome'

# Options to pass to vg mod for pruning phase. (if empty list, phase skipped)
# The primary path will always be added back onto the pruned grpah
prune-opts: ['-p', '-l', '16', '-S', '-e', '5']

# Options to pass to 2nd vg mod for pruning phase.  The input will
# will be piped in from the prune-opts command above.  Will be 
# skipped if empty list
prune-opts-2: ['-S', '-l', '32']

# Options to pass to vg kmers.
kmers-opts: ['-g', '-B', '-k', '16', '-H', '1000000000', '-T', '1000000001']

# Options to pass to vg gcsa indexing
gcsa-opts: ['-X', '3', '-Z', '3000']

########################
### vg_map Arguments ###

# Number of reads per chunk to use when splitting up fastq.  
# Each chunk will correspond to a vg map job
reads-per-chunk: 10000000

# Treat input fastq as paired-end interleaved
interleaved: False

# Core arguments for vg mapping (do not include file names or -t/--threads)
# Note -i/--interleaved will be ignored. use the --interleaved option 
# on the toil-vg command line instead
map-opts: ['-M2', '-W', '500', '-u', '0', '-U', '-O', '-S', '50', '-a']

# Type of vg index to use for mapping (either 'gcsa-kmer' or 'gcsa-mem')
index-mode: gcsa-mem

#########################
### vg_call Arguments ###
# Overlap option that is passed into make_chunks and call_chunk
overlap: 2000

# Chunk size
call-chunk-size: 10000000

# Context expansion used for graph chunking
chunk_context: 50

# Options to pass to chunk_gam. (do not include file names or -t/--threads)
filter-opts: ['-r', '0.9', '-fu', '-s', '1000', '-o', '0', '-q', '15']

# Options to pass to vg pileup. (do not include file names or -t/--threads)
pileup-opts: ['-q', '10', '-a']

# Options to pass to vg call. (do not include file names or -t/--threads)
call-opts: ['']

# Options to pass to vg genotype. (do not include file names or -t/--threads)
genotype-opts: ['']

# Use vg genotype instead of vg call
genotype: False

#########################
### vcfeval Arguments ###
# Options to pass to rgt vcfeval. (do not include filenaems or threads or BED)
vcfeval-opts: []

# BED region file for vcfeval
vcfeval-bed-regions:

""")

whole_genome_config = textwrap.dedent("""
# Toil VG Pipeline configuration file (created by toil-vg generate-config)
# This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
# Edit the values in the configuration file and then rerun the pipeline: "toil-vg run"
# 
# URLs can take the form: "/", "s3://"
# Local inputs follow the URL convention: "/full/path/to/input.txt"
# S3 URLs follow the convention: "s3://bucket/directory/file.txt"
#
# Comments (beginning with #) do not need to be removed. 
# Command-line options take priority over parameters in this file.  
######################################################################################################################

###########################################
### Toil resource tuning                ###

# These parameters must be adjusted based on data and cluster size
# when running on anything other than single-machine mode

# TODO: Reduce number of parameters here.  Seems fine grained, especially for disk/mem
# option to spin off config files for small/medium/large datasets?   

# The following parameters assign resources to small helper jobs that typically don't do 
    # do any computing outside of toil overhead.  Generally do not need to be changed. 
misc-cores: 1
misc-mem: '1G'
misc-disk: '1G'

# Resources allotted for xg indexing.
# this stage generally cannot take advantage of more than one thread
# For whole genome, suggest 200G mem and disk
xg-index-cores: 1
xg-index-mem: '200G'
xg-index-disk: '200G'

# Resources allotted for gcsa pruning.  Note that the vg mod commands used in
# this stage generally cannot take advantage of more than one thread
# For whole genome, suggest 60G mem and disk
prune-cores: 2
prune-mem: '60G'
prune-disk: '60G'

# Resources allotted for gcsa kmers.  
# For whole genome, suggest 70G mem and disk
kmers-cores: 16
kmers-mem: '70G'
kmers-disk: '60G'

# Resources allotted gcsa indexing
# For whole genome, suggest 200G mem and 3T disk
gcsa-index-cores: 32
gcsa-index-mem: '200G'
gcsa-index-disk: '3000G'

# Resources for fastq splitting and gam merging
# Important to assign as many cores as possible here for large fastq inputs
fq-split-cores: 32
fq-split-mem: '4G'
fq-split-disk: '200G'

# Number of threads to use for Rocksdb GAM indexing
# Generally, this should be kept low as speedup drops off radically 
# after a few threads.
gam-index-cores: 6

# Resources for *each* vg map job
# the number of vg map jobs is controlled by reads-per-chunk (below)
# For whole genome, suggest 100G memory and disk for whole genome
alignment-cores: 32
alignment-mem: '100G'
alignment-disk: '100G'

# Resources for chunking up a graph/gam for calling (and merging)
# typically take xg for whoe grpah, and gam for a chromosome,
# and split up into chunks of call-chunk-size (below)
# For whole genome, suggest 20G memory and 100G disk
call-chunk-cores: 2
call-chunk-mem: '20G'
call-chunk-disk: '100G'

# Resources for calling each chunk (currently includes pileup/call/genotype)
# For whole genome, suggest 20G memory and 20G disk
calling-cores: 32
calling-mem: '32G'
calling-disk: '32G'

# Resources for vcfeval
vcfeval-cores: 32
vcfeval-mem: '64G'
vcfeval-disk: '64G'

###########################################
### Arguments Shared Between Components ###
# Use output store instead of toil for all intermediate files (use only for debugging)
force-outstore: False

#############################
### Docker Tool Arguments ###
# Do not use docker for any commands
no-docker: False

## Docker Tool List ##
##   Each tool is specified as a list where the first element is the docker image URL,
##   and the second element indicates if the docker image has an entrypoint or not
##   If left blank or commented, then the tool will be run directly from the command line instead
##   of through docker. no-docker (above) overrides all these options. 

# Docker container to use for vg
vg-docker: ['quay.io/glennhickey/vg:v1.4.0-2213-gdce81f8', True]

# Docker container to use for bcftools
bcftools-docker: ['quay.io/cmarkello/bcftools', False]

# Docker container to use for tabix
tabix-docker: ['quay.io/cmarkello/htslib:latest', False]

# Docker container to use for jq
jq-docker: ['devorbitus/ubuntu-bash-jq-curl', False]

# Docker container to use for rtg
rtg-docker: ['realtimegenomics/rtg-tools:3.7.1', True]

# Docker container to use for pigz
pigz-docker: ['quay.io/glennhickey/pigz:latest', True]

##########################
### vg_index Arguments ###

# Name of index output files.  ex <name>.xg, <name>.gcsa etc. 
index-name: 'genome'

# Options to pass to vg mod for pruning phase. (if empty list, phase skipped)
# The primary path will always be added back onto the pruned grpah
prune-opts: ['-p', '-l', '16', '-S', '-e', '4']

# Options to pass to 2nd vg mod for pruning phase.  The input will
# will be piped in from the prune-opts command above.  Will be 
# skipped if empty list
prune-opts-2: ['-S', '-l', '32']

# Options to pass to vg kmers.
kmers-opts: ['-g', '-B', '-k', '16', '-H', '1000000000', '-T', '1000000001']

# Options to pass to vg gcsa indexing
gcsa-opts: ['-X', '3', '-Z', '3000']

########################
### vg_map Arguments ###

# Number of reads per chunk to use when splitting up fastq.  
# Each chunk will correspond to a vg map job
reads-per-chunk: 50000000

# Treat input fastq as paired-end interleaved
interleaved: False

# Core arguments for vg mapping (do not include file names or -t/--threads)
# Note -i/--interleaved will be ignored. use the --interleaved option 
# on the toil-vg command line instead
map-opts: []

# Type of vg index to use for mapping (either 'gcsa-kmer' or 'gcsa-mem')
index-mode: gcsa-mem

#########################
### vg_call Arguments ###
# Overlap option that is passed into make_chunks and call_chunk
overlap: 2000

# Chunk size
call-chunk-size: 8000000

# Context expansion used for graph chunking
chunk_context: 50


# Options to pass to chunk_gam. (do not include file names or -t/--threads)
filter-opts: ['-r', '0.9', '-fu', '-s', '1000', '-o', '0', '-q', '15', '-D', '20', '-C', '999']

# Options to pass to vg pileup. (do not include file names or -t/--threads)
pileup-opts: ['-q', '10', '-a']

# Options to pass to vg call. (do not include file names or -t/--threads)
call-opts: ['']

# Options to pass to vg genotype. (do not include file names or -t/--threads)
genotype-opts: ['']

# Use vg genotype instead of vg call
genotype: False

#########################
### vcfeval Arguments ###
# Options to pass to rgt vcfeval. (do not include filenaems or threads or BED)
vcfeval-opts: []

# BED region file for vcfeval
vcfeval-bed-regions:

""")

def generate_config(whole_genome = False):
    return whole_genome_config if whole_genome is True else default_config

def apply_config_file_args(args):
    """
    Merge args from the config file and the parser, giving priority to the parser.
    """

    # turn --*_opts from strings to lists to be consistent with config file
    for x_opts in ['map_opts', 'call_opts', 'filter_opts', 'genotype_opts']:
        if x_opts in args.__dict__.keys() and type(args.__dict__[x_opts]) is str:
            args.__dict__[x_opts] = args.__dict__[x_opts].split(' ')
            # get rid of any -t or --threads while we're at it
            for t in ['-t', '--threads']:
                if t in args.__dict__[x_opts]:
                    pos = args.__dict__[x_opts].index(t)
                    del args.__dict__[x_opts][pos:pos+2]

    # If no config file given, we generate a default one
    if args.config is None:
        config = generate_config()
    else:
        require(os.path.exists(args.config), 'Config, {}, not found. Please run '
            '"toil-vg generate-config > {}" to create.'.format(args.config, args.config))    
        with open(args.config) as conf:
            config = conf.read()
                
    # Parse config
    parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(config).iteritems()}
    options = argparse.Namespace(**parsed_config)

    # Add in options from the program arguments to the arguments in the config file
    #   program arguments that are also present in the config file will overwrite the
    #   arguments in the config file
    for args_key in args.__dict__:
        # Add in missing program arguments to config option list and
        # overwrite config options with corresponding options that are not None in program arguments
        if (args.__dict__[args_key]) or (args_key not in  options.__dict__.keys()):
            options.__dict__[args_key] = args.__dict__[args_key]
            
    return options

def config_subparser(parser):
    """
    Create a subparser for config.  Should pass in results of subparsers.add_parser()
    """

    parser.add_argument("--whole_genome", action="store_true",
        help="Make config tuned to process a whole genome on 32-core instances")


def config_main(args):
    """ config just prints out a file """
    
    sys.stdout.write(generate_config(args.whole_genome))
