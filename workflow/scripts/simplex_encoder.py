
from itertools import product
import numpy as np
import argparse
import sys

class WmerEncoder:

    def __init__(self):
        pass
    def encode(self):
        pass

class WMerSimplexEncoder(WmerEncoder):

    nuc_encodings = { 'A': [1,-1,-1], 'C': [-1,1,-1], 'G': [-1,-1,1], 'T':[1,1,1] }

    @classmethod
    def test(cls):
        test = cls(2)
        compare = np.array([[+1,-1,-1,-1,+1,+1,-1,+1,+1],
             [-1,+1,-1,+1,-1,+1,+1,-1,+1],
             [-1,-1,+1,+1,+1,-1,+1,+1,-1],
             [+1,+1,+1,-1,-1,-1,-1,-1,-1],
             [-1,+1,+1,+1,-1,-1,-1,+1,+1],
             [+1,-1,+1,-1,+1,-1,+1,-1,+1],
             [+1,+1,-1,-1,-1,+1,+1,+1,-1],
             [-1,-1,-1,+1,+1,+1,-1,-1,-1],
             [-1,+1,+1,-1,+1,+1,+1,-1,-1],
             [+1,-1,+1,+1,-1,+1,-1,+1,-1],
             [+1,+1,-1,+1,+1,-1,-1,-1,+1],
             [-1,-1,-1,-1,-1,-1,+1,+1,+1],
             [+1,-1,-1,+1,-1,-1,+1,-1,-1],
             [-1,+1,-1,-1,+1,-1,-1,+1,-1],
             [-1,-1,+1,-1,-1,+1,-1,-1,+1],
             [+1,+1,+1,+1,+1,+1,+1,+1,+1]])

        hopefully_zero = np.sum(test.oligonucleotide_encodings - compare)
        
        assert(hopefully_zero == 0)

    def __init__(self, w_mer_size = 3):

        #get all combinations (4^w) of nucleotides that may form a w-mer
        #for 2-mer: (A,A),(A,C) ... (T,G),(T,T)
        self.nucleotide_phrases = list(product(['A','C','G','T'], repeat = w_mer_size))

        nuc_by_position = list(zip(*self.nucleotide_phrases))
        #build monomer encoding matrix of shape (w, 4^w, 3), where the first index is the postion in the w-mer,
        #and the second two positions are the b1/b2 matrices for the paper
        self.monomer_encoding = np.array([
            [self.nuc_encodings[nuc] for nuc in nuc_by_position[w]]
            for w in range(w_mer_size)
        ])

        if w_mer_size > 1:

            #transpose so shape (w, 3, 4^w)
            self.monomer_encoding = np.transpose(self.monomer_encoding, (0,2,1))

            #get all possible combinations of outer products of different nucleotides
            self.product_indices = list(product(range(3), repeat = w_mer_size))

            #build outer product matrix from mononucleotide encodings
            #for w-mers, produces matrix shape (4^w, 3^w)
            self.oligonucleotide_encodings = []
            for outer_product_combo in self.product_indices: 
                #get list of accessor index pairs for (monomer nums, wky id)
                mononuc_index_combos = list(zip(*enumerate(outer_product_combo)))

                self.oligonucleotide_encodings.append(
                    np.multiply(*self.monomer_encoding[tuple(mononuc_index_combos)])
                )

            self.oligonucleotide_encodings = np.array(self.oligonucleotide_encodings).T

        else:
            self.oligonucleotide_encodings = np.squeeze(self.monomer_encoding)

        self.encoding_key = {''.join(nucleotides) : i for i, nucleotides in enumerate(self.nucleotide_phrases)}

    def encode(self, w_mers):
        try:
            return np.concatenate([
                self.oligonucleotide_encodings[self.encoding_key[w_mer], :] for w_mer in w_mers
            ])
        except KeyError:
            raise KeyError('W-mers {} could not be encoded, likely because of non-nucleotide character'.format(', '.join(w_mers)))


class WMerDummyEncoder(WmerEncoder):

    def __init__(self, w_mer_size = 2):
        
        self.nucleotide_phrases = list(product(['A','C','G','T'], repeat = w_mer_size))
        self.encoding_key = {''.join(nucleotides) : i for i, nucleotides in enumerate(self.nucleotide_phrases)}

    def encode(self, w_mers):
        
        dummy_encoding = []
        for w_mer in w_mers:
            try:
                dummy_vars = np.zeros(len(self.encoding_key))
                dummy_vars[self.encoding_key[w_mer]] = 1
                dummy_vars = dummy_vars[np.newaxis, :]
                dummy_encoding.append(dummy_vars)
            except KeyError:
                raise KeyError('W-mers {} could not be encoded, likely because of non-nucleotide character'.format(', '.join(w_mers)))
        
        return np.hstack(dummy_encoding)


class OligoEncoder:

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

    def __init__(self, oligo_length, wmer_encoder):
        self.encoders = [wmer_encoder(w) for w in range(1, oligo_length + 1)]
        self.oligo_length = oligo_length

    @staticmethod
    def find_ngrams(input_list, n):
        return zip(*[input_list[i:] for i in range(n)])


    def reverse_complement(self, seq):
        minus_strand = [self.complement[nuc] for nuc in seq]
        return ''.join(reversed(minus_strand))

    def encode(self, oligo):
        
        encoding = []

        for nuc_length, encoder in zip(range(1, self.oligo_length + 1), self.encoders):
            w_mers = [''.join(w_mer) for w_mer in self.find_ngrams(oligo, nuc_length)]
            encoding.append(encoder.encode(w_mers))

        return np.hstack(encoding)

    def get_nucleotide_components(self, oligo):
        components = []
        for i in range(1, self.oligo_length + 1):
            w_mers = [''.join(w_mer) for w_mer in self.find_ngrams(oligo, i)]
            components.extend(w_mers)
        return components

    def get_features(self, oligo):

        oligo = oligo.strip().upper()
        
        oligo_encoding = np.concatenate([
            self.encode(oligo)[np.newaxis, :], 
            self.encode(self.reverse_complement(oligo))[np.newaxis, :]
        ], axis = 1).astype(np.int8)

        return oligo_encoding
        

    def generate_feature_matrix(self, sequences):
        encodings = []
        oligo_length = None
        for oligo in sequences:
            if oligo_length is None:
                oligo_length = len(oligo)
            else:
                assert(len(oligo) == oligo_length)
            encodings.append(self.get_features(oligo))

        encodings = np.vstack(encodings)

        return encodings


if __name__ == "__main__":

    WMerSimplexEncoder.test()

    parser = argparse.ArgumentParser()
    parser.add_argument('oligos', nargs = "?", default = sys.stdin, type = argparse.FileType('r'))
    parser.add_argument('-o', '--output', required = True, type = str)
    parser.add_argument('-w','--wmer_len', type = int, default= 2)
    parser.add_argument('--header', action = 'store_true', default=False)
    parser.add_argument('--dummy',action='store_true',default=False)

    args = parser.parse_args()

    encoder = OligoEncoder(args.wmer_len, WMerSimplexEncoder if not args.dummy else WMerDummyEncoder)

    sequences = args.oligos.readlines()[1 if args.header else 0: ]

    encodings = encoder.generate_feature_matrix(sequences)

    print('Created {}x{} feature matrix!'.format(*[str(x) for x in encodings.shape]), file = sys.stderr)

    np.save(args.output, encodings)


