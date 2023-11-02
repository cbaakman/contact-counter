#!/usr/bin/env python

from argparse import ArgumentParser

import pandas
from matplotlib import pyplot


arg_parser = ArgumentParser()
arg_parser.add_argument("matrix_path")
arg_parser.add_argument("png_path")


if __name__ == "__main__":

    args = arg_parser.parse_args()

    data = pandas.read_csv(args.matrix_path)
    data = data.set_index("\\")

    figure = pyplot.figure(figsize=(13, 10))
    plot = figure.add_subplot()

    heatmap = plot.imshow(data, cmap="Greys", aspect="auto")
    figure.colorbar(heatmap)

    pyplot.xticks(range(len(data)), data.columns)
    pyplot.yticks(range(len(data)), data.index)

    figure.savefig(args.png_path, format="png")
