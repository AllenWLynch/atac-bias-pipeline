
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
            ops[op](group[colnum - 1])
        )

    return agg_results    

def iter_groups(*,fileinput, groupby, columns, ops, delim):

    last_index = None
    group = []
    for i, record in enumerate(fileinput):

        fields = record.split(delim)
        index = fields[groupby]

        if not index == last_index:

            if not index > last_index:
                raise NotSortedError('Record number {} in input file is out-of-order.'.format(str(i)))
            
            agg_results = release_group(group, columns, ops)
            yield (index, agg_results)

            last_index = index
            group = []

        group.append(fields)

    yield (index, release_group(group, columns, ops))

        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs = "?", default=sys.stdin)
    parser.add_argument('-g','--groupby',type = int, default=1)
    parser.add_argument('-o','--operations',default='sum',nargs='+', choices=['sum','mean','count','sample'])
    parser.add_argument('-c','--columns',type = int, default=2, nargs='+')
    parser.add_argument('--sample_num', type = int, default=150)
    parser.add_argument('-d','--delim',type = str, default='\t')
    args = parser.parse_args()

    print(args)
    assert(False)

    OPS = {
        'sum' : np.sum,
        'mean' : np.nanmean,
        'count' : np.count,
        'sample' : lambda x : np.random.choice(x, size = (min(len(x), args.sample_num),))
    }

    try:
        for index, results in iter_groups(fileinput=args.input, groupby=args.groupby, columns=args.columns,
                ops = [OPS[op] for op in args.ops], delim=args.delim):

            print(index,*results, sep = args.delim, file = sys.stderr)
    except NotSortedError as err:
        print(err, file = sys.stdout)
        exit(1)