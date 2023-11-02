#!/usr/bin/env python

import sys
import logging
import tarfile
from argparse import ArgumentParser
from typing import Dict, List
import subprocess
from threading import Thread
import csv
from math import log as ln
from tempfile import mkstemp

import torch
from Bio.PDB.Polypeptide import standard_aa_names

from counter.models.matrix import Matrix


_log = logging.getLogger(__name__)


arg_parser = ArgumentParser()
arg_parser.add_argument("exe")
arg_parser.add_argument("file_list", help="path to a file, containing paths to PDB files where contacts should be counted from")
arg_parser.add_argument("worker_count", type=int)
arg_parser.add_argument("output_file")


def count_for_one_structure(exe_path: str, pdb_path: str, chain1: str, chain2: str, cutoff_distance: float, device: torch.device) -> Matrix:

    result = subprocess.run([exe_path, pdb_path, str(cutoff_distance)], check=True, capture_output=True)

    m = Matrix(device=device)

    res_pairs = set([])
    for line in result.stdout.decode("utf_8").split("\n"):
        if len(line.strip()) == 0:
            continue

        # parse contact line
        aa1, _chain1, resnum1, name1, nr1, aa2, _chain2, resnum2, name2, nr2, distance = line.split()

        # filter chains involved
        if _chain1 != chain1 and _chain1 != chain2:
            continue
        if _chain2 != chain2 and _chain2 != chain1:
            continue

        # see if this residue pair is new
        res_pair = (chain1, int(resnum1), chain2, int(resnum2))
        if res_pair not in res_pairs:

            # if so, count it
            m.count_one(aa1, aa2)

        # remember which pairs we've already seen
        res_pairs.add(res_pair)

    return m


def extract_file(tar_path: str, file_path: str) -> str:
    tmp_file, tmp_path = mkstemp()
    os.close(tmp_file)

    with tarfile.open(tar_path, 'r') as tf:
        with tf.extractfile(file_path) as f:
            with open(tmp_path, 'wb') as output:
                output.write(f.read())

    return tmp_path


class CountThread(Thread):
    def __init__(
        self,
        exe_path: str,
        device: torch.device,
        file_list: List[str]
    ):
        Thread.__init__(self)

        self._device = device

        self.sum_of_matrices = Matrix(device=device)

        self._file_list = file_list

        self._exe_path = exe_path

    def run(self):
        for path in self._file_list:
            tmp_path = None
            try:
                tar_path, pdb_path = path.split(':')
                tmp_path = extract_file(tar_path, pdb_path)

                self.sum_of_matrices += count_for_one_structure(self._exe_path, tmp_path, "M", "P", 5.0, self._device)
            except:
                _log.exception(f"on {path}")
            finally:
                if tmp_path is not None:
                    os.remove(tmp_path)


if __name__ == "__main__":

    logging.basicConfig(stream=sys.stdout, level=logging.INFO)

    args = arg_parser.parse_args()

    device = torch.device("cpu")
    if torch.cuda.is_available():
        device = torch.device("cuda")

    # read pathnames from file
    with open(args.file_list, 'rt') as f:
        paths = f.read().strip().split('\n')

    # buld multithreads
    file_count = len(paths)
    file_frac = max(1, int(file_count / args.worker_count))

    threads = []
    file_index = 0
    for _ in range(args.worker_count):
        last_index = min(file_count, file_index + file_frac)

        t = CountThread(args.exe, device, paths[file_index: last_index])
        file_index = last_index

        t.start()
        threads.append(t)

    # count in multithreads
    sum_of_matrices = Matrix(device=device)
    for t in threads:
        t.join()
        sum_of_matrices += t.sum_of_matrices

    total = sum_of_matrices.sum()
    matrix_dict = sum_of_matrices.to_dict()

    # write matrix to output file
    with open(args.output_file, 'wt') as output_file:
        w = csv.writer(output_file)

        w.writerow(['\\'] + list(standard_aa_names))

        for aai in standard_aa_names:
            row = [aai]
            for aaj in standard_aa_names:
                p = (1 + matrix_dict[aai][aaj]) / total
                if p > 0.0:
                    e = -ln(p)
                else:
                    e = ""  # empty cell, NaN in pandas
                row.append(e)

            w.writerow(row)

