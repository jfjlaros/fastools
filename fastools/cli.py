import sys

from argparse import ArgumentParser, FileType, RawDescriptionHelpFormatter

from . import doc_split, version, usage
from .fastools import *
from .peeker import Peeker


def main():
    """Main entry point."""
    input_parser = ArgumentParser(add_help=False)
    input_parser.add_argument(
        'input_handle', metavar='INPUT', type=FileType('r'),
        help='input file')

    input2_parser = ArgumentParser(add_help=False)
    input2_parser.add_argument(
        'input_handles', metavar='INPUT', type=FileType('r'), nargs=2,
        help='input files')

    output_parser = ArgumentParser(add_help=False)
    output_parser.add_argument(
        'output_handle', metavar='OUTPUT', type=FileType('w'),
        help='output file')

    output2_parser = ArgumentParser(add_help=False)
    output2_parser.add_argument(
        'output_handles', metavar='OUTPUT', type=FileType('w'), nargs=2,
        help='output files')

    file_parser = ArgumentParser(
        add_help=False, parents=[input_parser, output_parser])

    qual_parser = ArgumentParser(add_help=False)
    qual_parser.add_argument(
        '-q', dest='quality', type=int, default=40,
        help='quality score (%(type)s default=%(default)s)')

    seq_parser = ArgumentParser(add_help=False)
    seq_parser.add_argument(
        'sequence', metavar='SEQ', type=str, help='a sequence (%(type)s)')

    range_parser = ArgumentParser(add_help=False)
    range_parser.add_argument(
        'first', metavar='FIRST', type=int,
        help='first base of the selection (%(type)s)')
    range_parser.add_argument(
        'last', metavar='LAST', type=int,
        help='last base of the selection (%(type)s)')

    name_parser = ArgumentParser(add_help=False)
    name_parser.add_argument(
        'name', metavar='ACCNO', type=str, help='accession number')

    description_parser = ArgumentParser(add_help=False)
    description_parser.add_argument(
        'description', metavar='DESCR', type=str,
        help='descriptino of the DNA sequence')

    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action='version', version=version(parser.prog))
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True

    subparser = subparsers.add_parser(
        'add', parents=[file_parser, seq_parser, qual_parser],
        description=doc_split(add))
    subparser.set_defaults(func=add)

    subparser = subparsers.add_parser(
        'aln', parents=[input2_parser], description=doc_split(aln))
    subparser.set_defaults(func=aln)

    subparser = subparsers.add_parser(
        'cat', parents=[input_parser], description=doc_split(cat))
    subparser.set_defaults(func=cat)

    subparser = subparsers.add_parser(
        'collapse', parents=[file_parser],
        description=doc_split(collapse))
    subparser.add_argument(
        '-s', '--stretch', dest='max_stretch', default=3, type=int,
        help='Length of the stretch (%(type)s default: %(default)s)')
    subparser.set_defaults(func=collapse)

    subparser = subparsers.add_parser('csv2fa2',
        parents=[input_parser, output2_parser], description=doc_split(csv2fa2))
    subparser.add_argument('-s', dest='skip_header', action='store_true',
        help='skip the first line of the CSV file')
    subparser.set_defaults(func=csv2fa2)

    subparser = subparsers.add_parser(
        'descr', parents=[input_parser], description=doc_split(descr))
    subparser.set_defaults(func=descr)

    subparser = subparsers.add_parser(
        'dna2rna', parents=[file_parser], description=doc_split(dna2rna))
    subparser.set_defaults(func=dna2rna)

    subparser = subparsers.add_parser(
        'edit', parents=[input_parser, output_parser],
        description=doc_split(edit))
    subparser.add_argument(
        'edits_handle', metavar='EDITS', type=FileType('r'),
        help='FASTA file containing edits')
    subparser.set_defaults(func=edit)

    subparser = subparsers.add_parser(
        'fa2fq', parents=[file_parser, qual_parser],
        description=doc_split(fa2fq))
    subparser.set_defaults(func=fa2fq)

    subparser = subparsers.add_parser(
        'fa2gb', parents=[file_parser, name_parser],
        description=doc_split(fa2gb))
    subparser.set_defaults(func=fa2gb)

    subparser = subparsers.add_parser(
        'famotif2bed', parents=[file_parser],
        description=doc_split(famotif2bed))
    subparser.add_argument(
        'motif', metavar='MOTIF', type=str, help='The sequence to be found')
    subparser.set_defaults(func=famotif2bed)

    subparser = subparsers.add_parser(
        'fq2fa', parents=[file_parser], description=doc_split(fq2fa))
    subparser.set_defaults(func=fq2fa)

    subparser = subparsers.add_parser(
        'gb2fa', parents=[file_parser], description=doc_split(gb2fa))
    subparser.set_defaults(func=gb2fa)

    subparser = subparsers.add_parser(
        'gen', parents=[output_parser, name_parser, description_parser],
        description=doc_split(gen))
    subparser.add_argument(
        'length', metavar='LENGTH', type=int,
        help='length of the DNA sequence')
    subparser.set_defaults(func=gen)

    subparser = subparsers.add_parser(
        'get', parents=[output_parser, name_parser],
        description=doc_split(get))
    subparser.add_argument(
        'email', metavar='EMAIL', type=str, help='email address')
    subparser.add_argument(
        '-s', dest='start', type=int, help='start of the area of interest')
    subparser.add_argument(
        '-p', dest='stop', type=int, help='end of the area of interest')
    subparser.add_argument(
        '-o', dest='orientation', type=int,
        help='orientation (1=forward, 2=reverse)')
    subparser.set_defaults(func=get)

    subparser = subparsers.add_parser(
        'length', parents=[input_parser], description=doc_split(length))
    subparser.set_defaults(func=length)

    subparser = subparsers.add_parser(
        'lenfilt', parents=[input_parser, output2_parser],
        description=doc_split(lenfilt))
    subparser.add_argument(
        '-l', dest='length', type=int, default=25,
        help='length threshold (%(type)s default: %(default)s)')
    subparser.set_defaults(func=lenfilt)

    subparser = subparsers.add_parser(
        'list_enzymes', description=doc_split(list_enzymes))
    subparser.set_defaults(func=list_enzymes)

    subparser = subparsers.add_parser(
        'maln', parents=[input_parser], description=doc_split(maln))
    subparser.set_defaults(func=maln)

    subparser = subparsers.add_parser(
        'mangle', parents=[file_parser], description=doc_split(mangle))
    subparser.set_defaults(func=mangle)

    subparser = subparsers.add_parser(
        'merge', parents=[input2_parser, output_parser],
        description=doc_split(merge))
    subparser.add_argument(
        '-f', dest='fill', type=int, default=0,
        help="Add 'N's between the reads (%(type)s default: %(default)s)")
    subparser.set_defaults(func=merge)

    subparser = subparsers.add_parser(
        'raw2fa',
        parents=[input_parser, output_parser, name_parser, description_parser],
        description=doc_split(raw2fa))
    subparser.set_defaults(func=raw2fa)

    subparser = subparsers.add_parser(
        'restrict', parents=[input_parser], description=doc_split(restrict))
    subparser.add_argument(
        '-r', dest='enzymes', metavar='ENZYME', type=str, action='append',
        default=[],
        help='restriction enzyme (use multiple times for more enzymes)')
    subparser.set_defaults(func=restrict)

    subparser = subparsers.add_parser(
        'reverse', parents=[file_parser], description=doc_split(reverse))
    subparser.set_defaults(func=reverse)

    subparser = subparsers.add_parser(
        'rna2dna', parents=[file_parser], description=doc_split(rna2dna))
    subparser.set_defaults(func=rna2dna)

    subparser = subparsers.add_parser(
        'rselect', parents=[file_parser, name_parser, range_parser],
        description=doc_split(rselect))
    subparser.set_defaults(func=rselect)

    subparser = subparsers.add_parser(
        's2i', parents=[file_parser], description=doc_split(s2i))
    subparser.set_defaults(func=s2i)

    subparser = subparsers.add_parser(
        'sanitise', parents=[file_parser], description=doc_split(sanitise))
    subparser.set_defaults(func=sanitise)

    subparser = subparsers.add_parser(
        'select', parents=[file_parser, range_parser],
        description=doc_split(select))
    subparser.set_defaults(func=select)

    subparser = subparsers.add_parser(
        'splitseq', parents=[input_parser, output2_parser, seq_parser],
        description=doc_split(splitseq))
    subparser.set_defaults(func=splitseq)

    subparser = subparsers.add_parser(
        'tagcount', parents=[input_parser, seq_parser],
        description=doc_split(tagcount))
    subparser.add_argument(
        '-m', dest='mismatches', type=int, default=2,
        help='amount of mismatches allowed (%(type)s default=%(default)s)')
    subparser.set_defaults(func=tagcount)

    sys.stdin = Peeker(sys.stdin)

    try:
        args = parser.parse_args()
    except IOError as error:
        parser.error(error)

    if args.subcommand == 'aln':
        for i in aln(args.input_handles):
            sys.stdout.write('{} {} {}'.format(*i))

    elif args.subcommand == 'length':
        sys.stdout.write(
            ' '.join(map(lambda x: str(x), length(args.input_handle))))

    elif args.subcommand == 'list_enzymes':
        sys.stdout.write('\n'.join(list_enzymes()))

    elif args.subcommand == 'restrict':
        sys.stdout.write(' '.join(
            map(lambda x: str(x), restrict(args.input_handle, args.enzymes))))

    elif args.subcommand == 'collapse':
        sys.stdout.write('Collapsed {} stretches longer than {}.'.format(
            collapse_fasta(
                args.input_handle, args.output_handle, args.max_stretch),
            args.max_stretch))

    elif args.subcommand == 's2i':
        sys.stdout.write('converted {} records'.format(s2i(
            args.input_handle, args.output_handle)))

    elif args.subcommand == 'tagcount':
        sys.stdout.write(
            count_tags(args.input_handle, args.sequence, args.mismatches))

    elif args.subcommand == 'cat':
        sys.stdout.write('\n'.join(cat(args.input_handle)))

    elif args.subcommand == 'descr':
        sys.stdout.write('\n'.join(descr(args.input_handle)))

    else:
        try:
            args.func(**{k: v for k, v in vars(args).items()
                if k not in ('func', 'subcommand')})
        except ValueError as error:
            parser.error(error)


if __name__ == '__main__':
    main()
