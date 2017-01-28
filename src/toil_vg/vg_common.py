#!/usr/bin/env python2.7
"""
Shared stuff between different modules in this package.  Some
may eventually move to or be replaced by stuff in toil-lib.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, timeit, errno
from uuid import uuid4
import pkg_resources, tempfile, datetime

from toil.common import Toil
from toil.job import Job
from toil_lib.toillib import *
from toil_lib.programs import docker_call

def add_docker_tool_parse_args(parser):
    """ centralize shared docker options and their defaults """
    parser.add_argument("--no_docker", action="store_true",
                        help="do not use docker for any commands")
    parser.add_argument("--vg_docker", type=str,
                        help="dockerfile to use for vg")
    parser.add_argument("--bcftools_docker", type=str,
                        help="dockerfile to use for bcftools")
    parser.add_argument("--tabix_docker", type=str,
                        help="dockerfile to use for tabix")
    parser.add_argument("--jq_docker", type=str,
                        help="dockerfile to use for jq")

def add_common_vg_parse_args(parser):
    """ centralize some shared io functions and their defaults """
    parser.add_argument('--config', default=None, type=str,
                        help='Config file.  Use toil-vg generate-config to see defaults/create new file')
    
    parser.add_argument("--checkpoint_all", action="store_true",
                        help="checkpoint all intermediate files to out store too (use only for debugging)")
                        
    parser.add_argument("--chroms", nargs='+',
                        help="Name(s) of reference path in graph(s) (separated by space).  If --graphs "
                        " specified, must be same length/order as --chroms")

    
def get_docker_tool_map(options):
    """ convenience function to parse the above _docker options into a dictionary """

    dmap = dict()
    if not options.no_docker:
        dmap["vg"] = options.vg_docker
        dmap["bcftools"] = options.bcftools_docker
        dmap["tabix"] = options.tabix_docker
        dmap["bgzip"] = options.tabix_docker
        dmap["jq"] = options.jq_docker

    # to do: could be a good place to do an existence check on these tools

    return dmap
        
class DockerRunner(object):
    """ Helper class to centralize docker calling.  So we can toggle Docker
on and off in just one place.  to do: Should go somewhere more central """
    def __init__(self, docker_tool_map = {}):
        # this maps a command to its full docker name
        # example:  docker_tool_map['vg'] = 'quay.io/ucsc_cgl/vg:latest'
        self.docker_tool_map = docker_tool_map

    def has_tool(self, tool):
        # return true if we have an image for this tool
        return tool in self.docker_tool_map

    def call(self, args, work_dir = '.' , outfile = None, errfile = None,
             check_output = False, inputs=[]):
        """ run a command.  decide to use docker based on whether
        its in the docker_tool_map.  args is either the usual argument list,
        or a list of lists (in the case of a chain of piped commands)  """
        # from here on, we assume our args is a list of lists
        if len(args) == 0 or len(args) > 0 and type(args[0]) is not list:
            args = [args]
        # convert everything to string
        for i in range(len(args)):
            args[i] = [str(x) for x in args[i]]
        if args[0][0] in self.docker_tool_map:
            return self.call_with_docker(args, work_dir, outfile, errfile, check_output, inputs)
        else:
            return self.call_directly(args, work_dir, outfile, errfile, check_output, inputs)
        
    def call_with_docker(self, args, work_dir, outfile, errfile, check_output, inputs): 
        """ Thin wrapper for docker_call that will use internal lookup to
        figure out the location of the docker file.  Only exposes docker_call
        parameters used so far.  expect args as list of lists.  if (toplevel)
        list has size > 1, then piping interface used """

        RealTimeLogger.get().info("Docker Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()
        
        if len(args) == 1:
            # just one command, use regular docker_call interface
            # where parameters is an argument list not including command
            tool = str(self.docker_tool_map[args[0][0]][0])
            tools = None
            if len(self.docker_tool_map[args[0][0]]) == 2:
                if self.docker_tool_map[args[0][0]][1]:
                    # not all tools have consistant entrypoints. vg and rtg have entrypoints
                    # but bcftools doesn't. This functionality requires the config file to
                    # operate
                    parameters = args[0][1:]
                else:
                    parameters = args[0]
            else:
                if args[0][0] in ["vg", "rtg"]:
                # todo:  not all tools have consistent entrypoints.  for instance vg and rtg
                # have entrypoints but bcftools doesn't.  hack here for now, but need to
                # have configurable entrypoints (or force image consistency) down the road. 
                    parameters = args[0][1:]
                else:
                    parameters = args[0]
        else:
            # there's a pipe.  we use the different piping interface that
            # takes in paramters as a list of single-string commands
            # that include arguments
            tool = None
            tools = str(self.docker_tool_map[args[0][0]][0])
            parameters = [" ".join(x) for x in args]

        ret = docker_call(tool=tool, tools=tools, parameters=parameters,
                           work_dir=work_dir, outfile = outfile,
                           errfile = errfile,
                           check_output = check_output,
                           inputs=inputs)

        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealTimeLogger.get().info("Successfully docker ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))

        return ret

    def call_directly(self, args, work_dir, outfile, errfile, check_output, inputs):
        """ Just run the command without docker """

        RealTimeLogger.get().info("Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()

        # this is all that docker_call does with the inputs parameter:
        for filename in inputs:
            assert(os.path.isfile(os.path.join(work_dir, filename)))

        procs = []
        for i in range(len(args)):
            stdin = procs[i-1].stdout if i > 0 else None
            if i == len(args) - 1 and outfile is not None:
                stdout = outfile
            else:
                stdout = subprocess.PIPE
            procs.append(subprocess.Popen(args[i], stdout=stdout, stderr=errfile,
                                          stdin=stdin, cwd=work_dir))
            
        for p in procs[:-1]:
            p.stdout.close()

        output, errors = procs[-1].communicate()
        for i, proc in enumerate(procs):
            sts = proc.wait()
            if sts != 0:            
                raise Exception("Command {} returned with non-zero exit status {}".format(
                    " ".join(args[i]), sts))

        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealTimeLogger.get().info("Successfully ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))            

        if check_output:
            return output

def get_files_by_file_size(dirname, reverse=False):
    """ Return list of file paths in directory sorted by file size """

    # Get list of files
    filepaths = []
    for basename in os.listdir(dirname):
        filename = os.path.join(dirname, basename)
        if os.path.isfile(filename):
            filepaths.append(filename)

    # Re-populate list with filename, size tuples
    for i in xrange(len(filepaths)):
        filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

    return filepaths

def clean_toil_path(path):
    """ Try to make input path into something toil friendly """
    # local path
    if ':' not in path:
        return 'file://' + os.path.abspath(path)
    else:
        return path

def init_out_store(options, command):
    """
    Write a little bit of logging to the output store.
    
    Rely on IOStore to create the store if it doesn't exist
    as well as to check its a valid location
     
    This is the only point in the code now where IOStore is used.
    Otherwise it's been replaced by Toil's export functions.  
    """
    f = tempfile.NamedTemporaryFile(delete=True)
    now = datetime.datetime.now()
    f.write('{}\ntoil-vg {} version {}\nOptions:'.format(now, command,
                    pkg_resources.get_distribution('toil-vg').version))
    for key,val in options.__dict__.items():
        f.write('{}: {}\n'.format(key, val))
    f.flush()
    IOStore.get(options.out_store).write_output_file(f.name, 'toil-vg-cmd.txt'.format(
        str(now).replace(':', '-').replace(' ', '-')))
    f.close()

    
def import_to_store(toil, options, path, copy_to_out_store = None,
                    out_store_key = None):
    """
    Imports a path into the File or IO store
    
    Returns the id. 

    (This is a relic of the old force_outstore logic.  Keep 
    it around in case we ever want to do something in all imports)
    """
    return toil.importFile(clean_toil_path(path))

def write_to_store(job, options, path, copy_to_out_store = None,
                   out_store_key = None):
    """
    Writes path into the File or IO store (from options.out_store)

    Returns the id in job's file store

    Abstract all store writing here, so we can switch on checkpointing
    centrally.  If checkpointing on, files are duplicated in the out store

    By default options.checkpoint_all is used to enable checkpointing. 
    This will be over-ridden by the copy_to_out_store parameter 
    if the latter is not None
    """
    file_id = job.fileStore.writeGlobalFile(path)
    
    if copy_to_out_store is True or (copy_to_out_store is None and options.checkpoint_all is True):
        out_store = IOStore.get(options.out_store)
        key = os.path.basename(path) if out_store_key is None else out_store_key
        exp_path = os.path.join(options.out_store, key)
        # copied from IOStore.get():
        if exp_path[0] in '/.':
            exp_path = 'file://' + exp_path
        else:
            store_type, store_arguments = exp_path.split(":", 1)
            if store_type == 'aws':
                region, bucket_name = store_arguments.split(":", 1)
                exp_path = 's3://' + bucket_name
            elif store_type == 'azure':
                account, container = store_arguments.split(":", 1)
                exp_path = 'wasb://' + container

        # Calling job.fileStore.exportFile() too soon after job.fileStore.writeGlobalFile
        # often causes the file to disappear or get truncated.
        time.sleep(10)
        RealTimeLogger.get().info("Exporting {} to output store: {}".format(str(file_id), exp_path))
        job.fileStore.exportFile(file_id, exp_path)

    return file_id 


def read_from_store(job, options, file_id, path = None, copy_to_out_store = None):
    """
    Reads id from the File store into path

    (No longer useful now that we never read from outstore. Will leave 
    for now in case we ever want to add something to all file reads)
    """
    return job.fileStore.readGlobalFile(file_id, path)

def write_dir_to_store(job, options, path, copy_to_out_store = None):
    """
    Need directory interface for rocksdb indexes.  Want to avoid overhead
    of tar-ing up as they may be big.  Write individual files instead, and 
    keep track of the names as well as ids (returns list of name/id pairs)
    """
    out_pairs = []
    file_list = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    for f in file_list:
        f_id = write_to_store(job, options, f, copy_to_out_store = copy_to_out_store,
                              out_store_key = path.replace('/', '_'))
        out_pairs.append(os.path.basename(f), f_id)
    return out_pairs

def read_dir_from_store(job, options, name_id_pairs, path = None, copy_to_out_store = None):
    """
    Need directory interface for rocksdb indexes.  Want to avoid overhead
    of tar-ing up as they may be big.  Takes as input list of filename/id pairs
    and reads these into the local directory given
    """
    if not os.path.isdir(path):
        os.mkdir(path)

    for name, key in name_id_pairs:
        read_from_store(job, options, key, os.path.join(path, name),
                        copy_to_out_store = copy_to_out_store)
    

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    
    From http://biopython.org/wiki/Split_large_file
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch
