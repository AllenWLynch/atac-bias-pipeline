
import numpy as np
from scipy import sparse
import argparse
import sys

def main(barcodes, peaks, fragments, corrected = True, deduplicate=False):

    barcodes = [barcode.strip() for barcode in barcodes]
    barcode_ids = {barcode : id_ for id_, barcode in enumerate(barcodes)}

    peaks = [peak.strip() for peak in peaks]

    peak_ids = {peak.replace('\t', '_') : id_ for id_, peak in enumerate(peaks)}

    row, column, data = [],[],[]
    
    dupgroups_encountered = {}
    for fragment in fragments:

        fragment = fragment.strip().split('\t')
        barcode = fragment[3]
        peak = '_'.join(fragment[10:13])
        dupgroup = fragment[4]

        add_count = not deduplicate
        if deduplicate and (dupgroup == '-1' or not dupgroup in dupgroups_encountered):
            dupgroups_encountered[dupgroup] = True
            add_count = True

        if add_count and peak in peak_ids and barcode in barcode_ids:
            row.append(barcode_ids[barcode])
            column.append(peak_ids[peak])
            if corrected:
                data.append(float(fragment[-1]))
            else:
                data.append(1.0)

    counts = sparse.coo_matrix(
        (data, (row, column)),
        shape = (len(barcodes), len(peaks))
    ).tocsc()

    return counts

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--barcodes', type = argparse.FileType('r'), required = True)
    parser.add_argument('-p','--peaks',type=argparse.FileType('r'), required=True)
    parser.add_argument('-f','--fragments', type = argparse.FileType('r'), required=True)
    parser.add_argument('--corrected',action='store_true', default=False)
    parser.add_argument('--dedup',action='store_true',default=False)
    parser.add_argument('-o','--output',type = str, required=True)

    args = parser.parse_args()

    count_matrix = main(args.barcodes, args.peaks, args.fragments,
        corrected=args.corrected, deduplicate=args.dedup)

    print('Created count matrix:\n' + repr(count_matrix), file = sys.stderr)
    
    sparse.save_npz(args.output, count_matrix)


