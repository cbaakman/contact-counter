#!/usr/bin/env python

import os
import re
import sys
import tarfile
from argparse import ArgumentParser

import pandas


arg_parser = ArgumentParser()
arg_parser.add_argument("list_csv", help="CSV file, containing all the IDs of the complexes (IEdb format)")
arg_parser.add_argument("models_dir", help="directory to search for model tars, may be in subdirectories")
arg_parser.add_argument("output_file")


seq_id_pattern = re.compile(r"% *SEQ ID: +(\d+\.\d+) *$")


def get_best_pandora_model(tar_path: str) -> str:

    id_ = os.path.splitext(os.path.basename(tar_path))[0]

    with tarfile.open(tar_path, 'r') as tf:
        filenames = tf.getnames()

        best_identity = 0.0
        best_model_name = ""
        for filename in filenames:
            if filename.startswith(f"{id_}/{id_}.") and filename.endswith(".pdb"):

                identity = 0.0
                with tf.extractfile(filename) as f:
                    for line in f:
                        line = line.decode("ascii")

                        match = seq_id_pattern.search(line)
                        if match is not None:
                            identity = float(match.group(1))

                if identity > best_identity:
                    best_identity = identity
                    best_model_name = os.path.basename(filename)

    return best_model_name


def list_files_under(path: str):
    if os.path.isdir(path):
        paths = []
        for dirname in os.listdir(path):
            paths += list_files_under(os.path.join(path, dirname))
        return paths
    else:
        return [path]


if __name__ == "__main__":

    args = arg_parser.parse_args()

    # read the complex ids from the CSV
    list_table = pandas.read_csv(args.list_csv)
    ids = list(list_table["ID"][:])

    # list all available tar files.
    tar_paths = list_files_under(args.models_dir)

    for id_ in ids:
        # look for the complex id among the tar files
        for tar_path in tar_paths:
            if os.path.basename(tar_path) == f"{id_}.tar":
                # pick the model in the tar
                model_filename = get_best_pandora_model(tar_path)

                # append name of pdb file and name of tar it's in
                with open(args.output_file, 'at') as f:
                    f.write(f"{tar_path}:{model_filename}\n")

                break
        else:
            raise FileNotFoundError(f"{args.models_dir}...{id_}.tar")
