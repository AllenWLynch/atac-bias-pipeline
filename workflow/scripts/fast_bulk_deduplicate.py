

import argparse
import subprocess
import sys
import os
from collections import defaultdict

class NotSortedError(Exception):
    pass

class FragmentWriter:

    def __init__(self, mark = False):
        self.mark = mark
        self.fragment_memory = defaultdict(dict)
        self.last_chrom = ''
        self.last_nuc = '-1'
        self.dup_groups_index = 0

    def add_fragment(self, chrom, start, end, barcode):

        if chrom < self.last_chrom:
            raise NotSortedError()
        elif chrom == self.last_chrom and int(start) < int(self.last_nuc):
            raise NotSortedError()
        
        if len(self.fragment_memory) == 0:
            self.last_chrom = chrom
            self.last_nuc = start

        if (chrom, start) == (self.last_chrom, self.last_nuc):
            location_id = (chrom, start, end)
            #if location / barcode combination not yet recorded, add to memory
            if not barcode in self.fragment_memory[location_id]:
                self.fragment_memory[location_id][barcode] = True
                #self.fragment_memory[location_id].append(barcode)
        else:
            raise ValueError('Cannot commit new fragment to memory without flushing previous fragments')
    
    def flush_memory(self):

        #for location_id, barcodes in sorted(self.fragment_memory.items(), key = lambda fragment : int(fragment[0][2])):
        for location_id, barcodes in self.fragment_memory.items():
            barcodes = list(barcodes.keys())
            if len(barcodes) > 1 and self.mark:
                for barcode in barcodes:
                    print('\t'.join([*location_id, barcode, str(self.dup_groups_index)]))
                self.dup_groups_index+=1
            else:
                print('\t'.join([*location_id, barcodes[0], '-1' if len(barcodes) == 1 else str(len(barcodes))]))

        del self.fragment_memory
        self.fragment_memory = defaultdict(dict)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Deduplicate a fragment file at the bulk level')
    parser.add_argument('input',nargs = '?', default=sys.stdin, type = argparse.FileType('r'), help = 'fragment file input, in format "chr\tstart\tend\tbarcode\t[pcr_duplicates]')
    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument('--mark', action='store_true', default=False, help = 'Mark duplicates in fragment file, but do not remove. If fragment is unique, mark with -1, else mark with an integer denoting its duplicate group membership.')
    group.add_argument('--assign', action = 'store_true', default=False, help = 'Assign a duplicate fragment to a cell.')
    parser.add_argument('-q','--quiet',action = 'store_true',default=False, help = 'Run quietly')

    
    args = parser.parse_args()

    fragment_writer = FragmentWriter(args.mark)

    for linenum, fragment in enumerate(args.input.readlines()):
        
        try:
            (chrom,start,end,barcode) = fragment.strip().split('\t')[:4]
        except ValueError:
            print('\nERROR: line #{} has incorrect number of fields. Make sure this file is tab-delineated.'.format(str(linenum + 1)), file = sys.stderr)
            exit(1)
        
        try:
            
            fragment_writer.add_fragment(chrom, start, end, barcode)
        
        except ValueError:
            fragment_writer.flush_memory()
            fragment_writer.add_fragment(chrom, start, end, barcode)
        except NotSortedError:
            print('\nERROR: line #{} is out of order. This will give an incorrect deduplication. Use cmd "sort -k1,1 -k2,2n [fragment_file]" before bulk deduplication.'.format(str(linenum+1)),
                file = sys.stderr)
            exit(1)

        if not args.quiet and linenum % 1000 == 0:
            print('\rProcessed {:.2f} million fragments ...'.format(linenum / 1000000), end = '', file = sys.stderr)

    fragment_writer.flush_memory()
    print('',file = sys.stderr)
    