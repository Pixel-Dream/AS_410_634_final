import sys

# Input
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        # The map for all sequences
        sequences = {}
        current_header = None
        current_seq = []

        for line in file:
            # Remove leading and trailing characters
            line = line.strip()
            # Encounter a new header
            if line.startswith('>'):
                # Record the previous header and sequence data
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                    # New sequence data
                    current_seq = []
                # New header
                current_header = line[1:]
            # Sequence data
            else:
                current_seq.append(line.upper().replace(' ', ''))

        # Record the last header and sequence data
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences
