
import argparse
import os
import sys
import joblib
import pyfaidx
from predict_cutsite_bias import process_sequences
from simplex_encoder import WMerSimplexEncoder, OligoEncoder
import multiprocessing
import numpy as np
from functools import partial
import warnings
from utils import Preprocessor

START_ADJUSTMENT = 8
END_ADJUSTMENT = 9

class RegionSequenceIterator:

    def __init__(self, chrom, start, end, fasta_file):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.faidx = pyfaidx.Fasta(fasta_file, default_seq = 'N', strict_bounds = False, read_ahead=end-start+100, sequence_always_upper=True)

    def __iter__(self):

        for cutsite_nuc, cutsite_start, cutsite_end in zip(
            range(self.start, self.end), 
            range(self.start - START_ADJUSTMENT, self.end - START_ADJUSTMENT), 
            range(self.start + END_ADJUSTMENT, self.end + END_ADJUSTMENT)
        ):
            cutsite_seq = self.faidx[self.chrom][cutsite_start:cutsite_end]

            if cutsite_seq.seq == '':
                yield (self.chrom, cutsite_nuc, cutsite_nuc+1), "N"*(START_ADJUSTMENT+END_ADJUSTMENT), "N"*(START_ADJUSTMENT+END_ADJUSTMENT)
            else:
                yield (self.chrom, cutsite_nuc, cutsite_nuc+1), cutsite_seq.seq, cutsite_seq.reverse.complement.seq

    def get_gc(self):
        return self.faidx[self.chrom][self.start:self.end].gc

class SeqGrouper:

    def __init__(self, seq_iterator, chunk_size):
        self.seqs = seq_iterator
        self.chunk_size = chunk_size

    def __iter__(self):
        buffer = []
        for seq in self.seqs:
            buffer.append(seq)
            if len(buffer) == self.chunk_size:
                yield list(zip(*buffer))
                buffer = []
        
        if len(buffer) > 0:
            yield list(zip(*buffer))

def summarize_bias(bias_scores):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        high_bias = np.quantile(bias_scores, 0.8)
        return np.nanmean(bias_scores[bias_scores > high_bias])

def summarize_window(window_cutsites, plus_sequences, minus_sequences, **kwargs):

        plus_biases, minus_biases = process_sequences(plus_sequences, **kwargs), process_sequences(minus_sequences, **kwargs) 
        
        plus_aggregate_score = summarize_bias(plus_biases)
        minus_aggregate_score = summarize_bias(minus_biases)

        window_chr, window_start, window_end = window_cutsites[0][0], window_cutsites[0][1], window_cutsites[-1][2]

        return (window_chr, window_start, window_end), plus_aggregate_score, minus_aggregate_score


def process_region(region,*,window_size, fasta_file, cutsite_kwargs, fragment_model, fragment_length):

    region_sequence_interface = RegionSequenceIterator(*region, fasta_file)

    results = []
    for window_cutsites, plus_sequences, minus_sequences in SeqGrouper(region_sequence_interface, window_size):

        if len(window_cutsites) > 20:
            window_region, plus_aggregate_score, minus_aggregate_score = summarize_window(window_cutsites, plus_sequences, minus_sequences, **cutsite_kwargs)
            
            gc_content = region_sequence_interface.get_gc()

            if gc_content == 0:
                gc_content = 0.5

            features = np.nan_to_num(np.array([fragment_length, gc_content, plus_aggregate_score, minus_aggregate_score]).reshape(1,-1), nan=0.5)
            
            try:
                duprate = fragment_model.predict(features)[0]
                results.append((window_region, duprate))
            except ValueError:
                pass

    return results

def preprocess_record(record):

    (chrom, start, end) = record.strip().split('\t')[:3]
    start, end = int(start), int(end)

    return chrom, start, end

    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('regions', type = argparse.FileType('r'), default=sys.stdin, nargs='?',help='Bedfile of regions to build background model')
    parser.add_argument('-m','--cutsite_model', type = str)
    parser.add_argument('-f','--fragment_model', type = str)
    parser.add_argument('-w','--windowsize', type = int, default=100)
    parser.add_argument('-g','--genome',type=str, help='Fasta file for indexing')
    parser.add_argument('-c','--cores',type = int, required=True)
    parser.set_defaults(wmer_len = 2, expected_length = 390, fragment_length = 50)

    args = parser.parse_args()

    fasta = pyfaidx.Fasta(args.genome, default_seq = 'N', strict_bounds = False, read_ahead=100000)
    encoder = OligoEncoder(args.wmer_len, WMerSimplexEncoder)
    model = joblib.load(args.cutsite_model)
    fragment_model = joblib.load(args.fragment_model)

    process_partial = partial(process_region, window_size = args.windowsize, fasta_file = args.genome, 
            cutsite_kwargs = dict(model = model, encoder = encoder, expected_length = args.expected_length),
            fragment_model = fragment_model, fragment_length=args.fragment_length)

    windows_written = 0
    with multiprocessing.Pool(args.cores) as p:

        for region_results in p.imap(process_partial, map(preprocess_record, args.regions.readlines())):

            for region, duprate in region_results:

                print(*region, duprate, sep='\t',file=sys.stdout)
                windows_written+=1
                if windows_written % 50 == 0:
                    print('\rWindows written: {}'.format(str(windows_written)), end = '', file=  sys.stderr)
    
    print('',file = sys.stderr)
