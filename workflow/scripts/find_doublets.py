
import scrublet as scr
import argparse
from scipy import sparse
import sys
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('count_matrix', type = str)
    parser.add_argument('output', type = argparse.FileType('w'))
    args = parser.parse_args()

    count_matrix = sparse.load_npz(args.count_matrix)

    scrub = scr.Scrublet(count_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    for score, prediction in zip(doublet_scores, predicted_doublets):
        print(score, prediction, sep = '\t', file = args.output)

