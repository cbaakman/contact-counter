#!/usr/bin/env python

import os
import re
import sys
import tarfile
from argparse import ArgumentParser
from glob import glob

import pandas


arg_parser = ArgumentParser()
arg_parser.add_argument("list_csv", help="CSV file, containing all the IDs of the complexes (IEdb format)")
arg_parser.add_argument("models_dir", help="directory to search for model tars, may be in subdirectories")
arg_parser.add_argument("output_dir")


seq_id_pattern = re.compile(r"% *SEQ ID: +(\d+\.\d+) *$")


def get_best_pandora_model(tar_path: str) -> str:

    id_ = os.path.splitext(os.path.basename(tar_path))[0]

    with tarfile.open(tar_path, 'r') as tf:
        filenames = tf.getnames()

        best_identity = 0.0
        best_model = ""
        best_model_name = ""
        for filename in filenames:
            if filename.startswith(f"{id_}/{id_}.") and filename.endswith(".pdb"):

                identity = 0.0
                model_s = ""
                with tf.extractfile(filename) as f:
                    for line in f:
                        line = line.decode("ascii")
                        model_s += line

                        match = seq_id_pattern.search(line)
                        if match is not None:
                            identity = float(match.group(1))

                if identity > best_identity:
                    best_identity = identity
                    best_model = model_s
                    best_model_name = os.path.basename(filename)

    return best_model, best_model_name


def list_files_under(path: str):

    if os.path.isdir(path):
        for dirname in os.listdir(path):
            return list_files_under(os.path.join(path, dirname))

    else:
        return [path]


if __name__ == "__main__":

    args = arg_parser.parse_args()

    if os.path.isdir(args.output_dir):
        if len(os.listdir(args.output_dir)) > 0:
            raise ValueError(f"{args.output_dir} exists and isn't empty")
    else:
        os.mkdir(args.output_dir)

    list_table = pandas.read_csv(args.list_csv)
    ids = list(list_table["ID"][:])

    tar_paths = list_files_under(args.models_dir)

    paths = []
    for id_ in ids:
        for tar_path in tar_paths:
            if os.path.basename(tar_path) == f"{id_}.tar":
                model_s, model_filename = get_best_pandora_model(tar_path)
                break
        else:
            raise FileNotFoundError(f"{id_}.tar")

        output_path = os.path.join(args.output_dir, model_filename)
        with open(output_path, 'wt') as f:
            f.write(model_s)

        list_path = os.path.join(args.output_dir, "model_list.txt")
        with open(list_path, 'at') as f:
            f.write(output_path + "\n")
