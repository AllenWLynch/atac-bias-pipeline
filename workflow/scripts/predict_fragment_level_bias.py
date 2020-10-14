
from joblib import load
import numpy as np
import argparse
import sys

from utils import lines_grouper, Preprocessor


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bias_file',type = argparse.FileType('r'), default = sys.stdin, nargs = '?')
    parser.add_argument('-m','--fragment_model', type=str, required=True)
    parser.add_argument('--bias1',type=int, default = 6)
    parser.add_argument('--bias2', type = int, default=7)
    parser.add_argument('--gc', type = int, default=8)
    parser.add_argument('--len', type = int, default=9)
    parser.add_argument('--fillnan', action = 'store_true', help = 'Fill nan cutsite biases with previously-encountered value (in sorted fragment files, the biases of nearby sites may be similar)')

    args = parser.parse_args()

    print('Predicting bias ...', file = sys.stderr)
    
    pipeline = load(args.fragment_model)

    grouper = lines_grouper(args.bias_file, 1000)

    fragments_processed = 0
    for fragments in iter(grouper):

        fragments_processed += len(fragments)
        print('\rFragments processed: {}'.format(str(fragments_processed)), end = '', file = sys.stderr)
        
        fragment_features = list(zip(*[x.strip().split('\t') for x in fragments]))

        feature_matrix = np.vstack([
            fragment_features[args.len - 1],
            fragment_features[args.gc - 1],
            fragment_features[args.bias1 - 1],
            fragment_features[args.bias2 - 1],
        ]).astype(np.float32).T

        if np.isnan(feature_matrix).any():
            if args.fillnan:
                feature_matrix = np.where(np.isnan(feature_matrix), np.nanmean(feature_matrix, axis = 0), feature_matrix)
            else:
                raise AssertionError('Encountered nan value in features. Try using --fillnan to replace nan with mean baises.')

        try:
            dup_rates = pipeline.predict(feature_matrix)
        except ValueError:
            raise ValueError("Encountered nan value in features. Try using --fillnan to replace nan with mean baises. If this option is flagged, then your entire sample may contain nans for some feature.")


        for fragment, dup_rate in zip(fragments, dup_rates):
            print(fragment.strip(), str(dup_rate), sep = '\t')        
    
    print('', file = sys.stderr)