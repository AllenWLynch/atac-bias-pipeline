
import argparse
import os
import sys
import numpy as np
from functools import partial

class NotSortedError(Exception):
    pass

def release_group(group, columns, ops):

    group = list(zip(*group))

    agg_results = []
    for colnum, op in zip(columns, ops):

        agg_results.append(
            op(np.array(group[colnum - 1]).astype(np.float))
        )

    return agg_results    

def iter_groups(*,fileinput, groupby, columns, ops, delim, ignore_warnings = False):

    last_index = None
    group = []
    for i, record in enumerate(fileinput):

        fields = record.split(delim)
        index = fields[groupby - 1]

        if not last_index is None and not index == last_index:

            if ignore_warnings==False and index < last_index:
                raise NotSortedError('Record number {} in input file is out-of-order.'.format(str(i)))
            
            agg_results = release_group(group, columns, ops)
            yield (last_index, agg_results)
            group = []

        last_index = index
        group.append(fields)

    yield (last_index, release_group(group, columns, ops))

        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs = "?", type = argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-g','--groupby',type = int, default=1)
    parser.add_argument('-o','--operations',default='sum',nargs='+', choices=['sum','mean','count','sample','median'])
    parser.add_argument('-c','--columns',type = int, default=2, nargs='+')
    parser.add_argument('--sample_num', type = int, default=150)
    parser.add_argument('-d','--delim',type = str, default='\t')
    parser.add_argument('--ignore_warnings',action='store_true')
    args = parser.parse_args()

    OPS = {
        'sum' : np.sum,
        'mean' : np.nanmean,
        'count' : np.size,
        'sample' : lambda x : list(np.random.choice(x, size = (min(len(x), args.sample_num),), replace = False)),
        'median' : np.median,
    }

    try:
        for index, results in iter_groups(fileinput=args.input, groupby=args.groupby, columns=args.columns,
                ops = [OPS[op] for op in args.operations], delim=args.delim, ignore_warnings=args.ignore_warnings):

            print(index,*results, sep = args.delim, file = sys.stdout)

    except NotSortedError as err:
        print(err, file = sys.stderr)
        exit(1)