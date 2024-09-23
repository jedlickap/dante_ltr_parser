from Bio import SeqIO

class FastaParser:
    def __init__(self, file_path):
        self.file_path = file_path

    def get_longest_sequence_in_mb(self):
        longest_length = 0

        # Parse the multi-FASTA file
        for record in SeqIO.parse(self.file_path, "fasta"):
            seq_length = len(record.seq)
            if seq_length > longest_length:
                longest_length = seq_length

        # Convert length to megabase pairs (Mb)
        longest_length_mb = longest_length / 1_000_000

        return longest_length_mb

# Example usage:
# parser = FastaParser("path_to_your_file.fasta")
# longest_sequence, length_in_mb = parser.get_longest_sequence_in_mb()
# print(f"Longest sequence (in Mb): {length_in_mb}")
# print(f"Longest sequence: {longest_sequence}")
