import so3g
from spt3g import core
import numpy as np
import os
import sys
import csv
import argparse


_UNITS = {
    'bytes': 1,
    'kb': 1024,
    'mb': 1024*1024,
    'gb': 1024*1024*1024,
}

def get_parser():
    parser = argparse.ArgumentParser(
        epilog="Run '%(prog)s COMMAND --help' to see additional "
        "details and options.")
    cmdsubp = parser.add_subparsers(
        dest='mode')

    # Shared arguments for subprocessors ...
    data_args = argparse.ArgumentParser(add_help=False)
    data_args.add_argument(
        'files', nargs='+', default=None,
        help="One or more G3 files to scan.")
    data_args.add_argument(
        '--recursive', '-r', action='store_true',
        help="All arguments are traversed recursively; only files "
        ".g3 extension files are scanned.")
    data_args.add_argument(
        '--strip-tokens', default='observatory.feeds',
        help="Tokens to hide in provider and field names. "
        "Pass this is a single .-delimited string.")
    data_args.add_argument(
        '--block-size', '-B', default='b',
        help="Summarize storage in units of bytes,kB,MB,GB (pass b,k,M,G).")
    data_args.add_argument(
        '--sort-size', '-s', action='store_true',
        help="Sort results, if applicable, by size (descending).")

    output_args = argparse.ArgumentParser(add_help=False)
    output_args.add_argument(
        '--csv', help=
        "Store data as CSV to specified filename.")

    # Main "mode" subprocessors.

    # "list-files"
    p = cmdsubp.add_parser(
        'list-files',
        parents=[data_args, output_args],
        help="Report per-file stats.",
        usage="""

    %(prog)s [options] FILE [...]

        This module reads each file and reports basic stats such as size and
        whether the stream is valid.
        """)

    # "list-provs"
    p = cmdsubp.add_parser(
        'list-provs',
        parents=[data_args, output_args],
        help="List all data providers (feeds).",
        usage="""

    %(prog)s [options] FILE [...]

        This module reads all specified files and reports a list of
        all data providers (a.k.a. feeds) encountered in the data,
        along with total data volume and average frame size, per
        provider.
        """)

    # "list-fields"
    p = cmdsubp.add_parser(
        'list-fields',
        parents=[data_args, output_args],
        help="List all data field names.",
        usage="""

    %(prog)s [options] FILE [...]

        This module reads all specified files and reports a list of
        all data fields with their total sample count.
        """)

    # Done.
    return parser


def get_file_list(args, suffix='.g3'):
    if args.recursive:
        all_files = []
        for root in args.files:
            these_files = []
            for base, dirs, files in os.walk(root):
                for f in files:
                    if f.endswith(suffix):
                        these_files.append(os.path.join(base, f))
            all_files.extend(sorted(these_files))
        return all_files
    else:
        return args.files


def format_table(rows, header=None, fmts=None, align=None):
    """Return a string with data from rows organized into a text table.

    If you pass header, it must be a list with the same number of
    elements as the first row.

    If you pass fmts or align, it must be a dict where the key is
    either the column index or the corresponding header name.

    For fmts, each dict value must be an anonymous python format
    string, e.g. "{:5.1f}".  For align, each dict value must be either
    'right' or 'left'.

    Default alignment is based on data types in the first row of data
    -- if float or int, it will be right aligned.  Otherwise, left.

    """
    if header is None:
        if len(rows):
            ncol = len(rows[0])
    else:
        ncol = len(header)

    def dict_to_per_col(data, default=None, col_defaults=None):
        # Change dict into a list, per-column, taking default from
        # col_defaults dict if possible and then from default.
        if data is None:
            data = {}
        if col_defaults is None:
            col_defaults = [default] * ncol
        elif isinstance(col_defaults, dict):
            col_defaults = dict_to_per_col(col_defaults, default)
        output = []
        for i in range(ncol):
            if i in data:
                output.append(data[i])
            elif header is not None and header[i] in data:
                output.append(data[header[i]])
            else:
                output.append(col_defaults[i])
        return output

    fmts = dict_to_per_col(fmts, '{}')
    fmt_align = {}
    for i, fmt in enumerate(fmts):
        if '>' in fmt:
            fmt_align[i] = 'right'
        elif '<' in fmt:
            fmt_align[i] = 'left'
        elif len(rows) and isinstance(rows[0][i], (float, int)):
            fmt_align[i] = 'right'
        else:
            fmt_align[i] = 'left'
    align = dict_to_per_col(align, 'left', fmt_align)

    # First round formatting ...
    rows = [[fmt.format(c) for fmt, c in zip(fmts, row)]
            for row in rows]

    # Now set column widths ...
    col_lens = [max([0 if header is None else len(header[i])]
                    + [len(r[i]) for r in rows])
                for i in range(ncol)]

    # Second round formatting ...
    align_fmts = [('{:<%is}' if al == 'left' else '{:>%is}') % cl
             for al, cl in zip(align, col_lens)]
    rows = [[a.format(c) for a, c in zip(align_fmts, r)]
            for r in [header] + rows]

    # Insert hline?
    if header is not None:
        line = ['-' * len(c) for c in rows[0]]
        rows.insert(1, line)

    return '\n'.join([' '.join(r) for r in rows])


def convert_units(rows, cols, args):
    units = args.block_size.lower()
    if units == 'b':
        units = 'bytes'
    if units not in _UNITS:
        units += 'b'
    if units not in _UNITS:
        raise ValueError(f'Cannot interpret {args.block_size} '
                         'as a block size (b,k,M,G).')
    factor = _UNITS[units]
    if factor == 1:
        fmt = '%d'
    else:
        fmt = '%.1f'
    rows = [
        tuple([(fmt % (c/factor)) if i in cols else c
               for i, c in enumerate(r)])
        for r in rows]
    return rows, units


def write_csv(filename, rows, header=None):
    with open(filename, 'w', newline='') as csvfile:
        cw = csv.writer(csvfile)#, dialect='memberdb')
        if header is not None:
            cw.writerow(header)
        for r in rows:
            cw.writerow(list(r))

def produce_output(rows, header, fmts=None, align=None, csv=False):
    if csv:
        write_csv(csv, rows, header=header)
    else:
        print(format_table(rows, fmts=fmts, header=header, align=align))


class TokenCleanser:
    def __init__(self, boring_tokens):
        self.boring = [t for t in boring_tokens.split('.')
                       if len(t)]
    def __call__(self, name):
        if len(self.boring):
            name = '.'.join([t for t in name.split('.')
                             if t not in self.boring])
        return name


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = get_parser()
    args = parser.parse_args(args=args)

    if args.mode == 'list-files':
        rows = []
        file_list = get_file_list
        for filename in get_file_list(args):
            file_size = os.path.getsize(filename)
            clean_exit = True
            r = core.G3Reader(filename)
            while True:
                try:
                    f = r.Process(None)
                except:
                    clean_exit = False
                    break
                end = r.tell()
                if f is None or len(f) == 0:
                    break
            rows.append((filename, file_size, end,
                         {True: 'no', False: 'YES'}[clean_exit]))

        if args.sort_size:
            rows = sorted(rows, key=lambda row: -row[1])

        rows, units = convert_units(rows, [1, 2], args)

        header = ['filename', f'size_{units}',
                  f'usable_{units}', 'error']
        produce_output(rows, header, align={1: 'right', 2: 'right'},
                       csv=args.csv)

    elif args.mode == 'list-provs':
        counts = {}
        renamer = TokenCleanser(args.strip_tokens)

        for filename in get_file_list(args):
            r = core.G3Reader(filename)
            while True:
                start = r.tell()
                f = r.Process(None)
                end = r.tell()
                if f is None or len(f) == 0:
                    break
                f = f[0]
                if 'address' in f:
                    key = renamer(f['address'])
                    if key not in counts:
                        counts[key] = []
                    counts[key].append(end - start)

        rows = [(k, np.sum(v).tolist(), np.mean(v))
                for k, v in sorted(counts.items())]

        if args.sort_size:
            rows = sorted(rows, key=lambda row: -row[1])

        rows, units = convert_units(rows, [1, 2], args)
        header = ['provider_name', f'total_{units}', f'frame_{units}']
        produce_output(rows, header, align={1: 'right', 2: 'right'},
                       csv=args.csv)

    elif args.mode == 'list-fields':
        # Count samples.
        counts = {}
        renamer = TokenCleanser(args.strip_tokens)

        for filename in get_file_list(args):
            r = core.G3Reader(filename)
            while True:
                f = r.Process(None)
                if f is None or len(f) == 0:
                    break
                f = f[0]
                if 'address' in f:
                    addr = renamer(f['address'])
                    for block in f['blocks']:
                        keys = block.keys()
                        n = len(block.times)
                        for k in keys:
                            field = addr + '.' + k
                            if not field in counts:
                                counts[field] = 0
                            counts[field] += n
        header = ['field_name', 'samples']
        rows = sorted(counts.items())

        if args.sort_size:
            rows = sorted(rows, key=lambda row: -row[1])

        produce_output(rows, header, csv=args.csv)

    elif args.mode is None:
        parser.error('Provide a valid mode. (See --help for more.)')

    else:
        print(f'Unimplemented mode "{args.mode}"!')


if __name__ == '__main__':
    main()
