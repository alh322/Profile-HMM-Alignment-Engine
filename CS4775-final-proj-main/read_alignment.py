import argparse


def read_alignment(alignment_file, format):
    """
    Reads a file with specified format containing an alignment and returns a list of sequences
    """
    sequences = []
    if format == 'fasta':
        with open(alignment_file, 'r') as f:
            for line in f:
                if line[0] == '>':
                    sequences.append('')
                else:
                    sequences[-1] += line.strip()
    if format == 'msf':
        with open(alignment_file, 'r') as file:
            lines = file.readlines()

        start_processing = False
        sequences_temp = {}
        for line in lines:
            if '//' in line:
                start_processing = True
                continue
            if start_processing:
                parts = line.split()
                if parts:
                    seq_id = parts[0]
                    seq_data = ''.join(parts[1:])
                    sequences_temp[seq_id] = sequences_temp.get(
                        seq_id, '') + seq_data

        sequences = [s. replace('.', '-')
                     for s in list(sequences_temp.values())]
    return sequences


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads a file with an alignment and returns a list of sequences')
    parser.add_argument('alignment_file', type=str,
                        help='Path to the alignment file')
    parser.add_argument('format', type=str,
                        help='Format of the alignment file (fasta or align)')
    args = parser.parse_args()
    sequences = read_alignment(args.alignment_file, args.format)
    print(sequences)
