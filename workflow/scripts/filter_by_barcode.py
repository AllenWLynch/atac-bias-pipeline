
import argparse
import sys
import os

def get_allowed_barcodes(barcode_input):
    return {
        barcode : True for barcode in map(lambda x : x.strip(), barcode_input)
    }

def apply_filter(fragments, allowed_values, colnum):

    values_dict = get_allowed_barcodes(allowed_values)

    for fragment in fragments:
        fields = fragment.strip().split('\t')
        if fields[colnum] in values_dict:
            yield fragment.strip()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('fragments',type=argparse.FileType('r'))
    parser.add_argument('barcodes', type = argparse.FileType('r'))
    args = parser.parse_args()

    for fragment in apply_filter(args.fragments, args.barcodes, 3):
        print(fragment, file = sys.stdout)

