import argparse
import sys
import os

def dir_path(dir_path):
    if os.path.isdir(dir_path):
        return dir_path
    else:
        raise argparse.ArgumentTypeError(dir_path + ' is not a valid directory.')

def main(bedfile, dir_path, pattern):

    assert(r'{chrom}' in pattern), r'Pattern provided must contain wildcard "{chrom}"'

    chrom_files = {}
    try:
        for fragment in bedfile:

            chrom = fragment.split('\t')[0]

            if not chrom in chrom_files:
                chrom_files[chrom] = open(os.path.join(dir_path, pattern.format(chrom = chrom)), 'w')

            chrom_files[chrom].write(fragment)

    finally:
        for filehandle in chrom_files.values():
            filehandle.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('bedfile',type = argparse.FileType('r'), nargs= '?', default=sys.stdin)
    parser.add_argument('-d','--directory', type = dir_path, default = './')
    parser.add_argument('--pattern', type = str, required=False)
    args = parser.parse_args()

    pattern = args.pattern or '{chrom}_' + os.path.basename(args.bedfile.name)

    try:
        main(args.bedfile, args.directory, pattern)
    except Exception as err:
        print('ERROR:', err, file = sys.stderr)