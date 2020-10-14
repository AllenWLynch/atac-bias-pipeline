
from joblib import load
import numpy as np
import argparse
import sys
import warnings

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
    feature_means = None
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

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            if np.isnan(feature_matrix).any():
                if args.fillnan:
                    if feature_means is None:
                        feature_means = np.nanmean(feature_matrix, axis = 0)
                    else:
                        #replace feature_means with local fragment means
                        new_mean = np.nanmean(feature_matrix, axis = 0)
                        #replace previous mean with local mean if all columns have atleast one representative
                        if not np.isnan(new_mean).any():
                            feature_means = new_mean
                        feature_matrix = np.where(np.isnan(feature_matrix), feature_means, feature_matrix)
                else:
                    raise AssertionError('Encountered nan value in features. Try using --fillnan to replace nan with mean baises.')

        if not np.isnan(feature_matrix).any():
            dup_rates = pipeline.predict(feature_matrix)
        else:
            print('\n(non-fatal) ERROR: Encountered group of fragments with nan features that could not be filled with a group mean. Assigning duplication rate of "nan"', file = sys.stderr)
            dup_rates = np.full(feature_matrix.shape[0], np.nan)

        for fragment, dup_rate in zip(fragments, dup_rates):
            print(fragment.strip(), str(dup_rate), sep = '\t')        
    
    print('', file = sys.stderr)