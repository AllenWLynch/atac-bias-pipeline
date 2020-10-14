
import joblib
import numpy as np
import argparse
import sys
from simplex_encoder import OligoEncoder, WMerSimplexEncoder
import multiprocessing
from functools import partial

class lines_grouper:

    def __init__(self, file_object, chunk_size):
        self.file_object = file_object
        self.chunk_size = chunk_size

    def __iter__(self):
        morelines = True
        while morelines:
            group = []
            for i in range(self.chunk_size):
                nextline = self.file_object.readline()
                if nextline == '':
                    morelines = False
                    break
                else:
                    group.append(nextline)
            if len(group) > 0:
                yield group

def process_sequences(sequences, *, encoder, model, expected_length):

    #sequences = [x.strip() for x in sequences]
    features_list = []
    invalid_sequences = []
    for i, oligo in enumerate(sequences):
        try:
            features = encoder.get_features(oligo)

            if features.shape[-1] != expected_length:
                raise KeyError()

        except KeyError:
            #make up dummy features
            features = np.ones((1, expected_length))
            invalid_sequences.append(i)

        features_list.append(features)

    features_matrix = np.vstack(features_list)

    bias_predictions = model.predict_proba(features_matrix)[:,1]

    bias_predictions[invalid_sequences] = np.nan

    return bias_predictions


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--sequences',type = argparse.FileType('r'), default = sys.stdin, nargs = '?')
    parser.add_argument('-m','--model', type=str, required=True)
    parser.add_argument('-w','--wmer_len', type = int, default= 2)
    parser.add_argument('-c','--cores', type = int, default = 1)
    parser.add_argument('-e','--expected_length', type = int, required=True, help = 'Expected length of featurization or nucleotide sequences. Used to reject invalid nucleotide sequences.')

    args = parser.parse_args()

    print('Predicting bias ...', file = sys.stderr)
    encoder = OligoEncoder(args.wmer_len, WMerSimplexEncoder)

    model = joblib.load(args.model)

    sequences_processed = 0

    process_partial = partial(process_sequences, encoder = encoder, model = model, expected_length = args.expected_length)

    with multiprocessing.Pool(args.cores) as pool:

        for i, biases in enumerate(pool.imap(process_partial, lines_grouper(args.sequences, 10000))):
            
            sequences_processed+=len(biases)
            print('\rProcessed sequences: ', sequences_processed, file = sys.stderr, end = '')

            for bias in biases:
                print(str(bias), file = sys.stdout)

    print('', file = sys.stderr)