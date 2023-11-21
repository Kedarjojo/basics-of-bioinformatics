def revcomp(seq):

    comp_base={'A':'T','G':'C','T':'A','C':'G'}
    rev_seq = seq[::-1]
    rev_comp = ''

    for base in rev_seq:
        rev_comp+=comp_base[base]
    return rev_comp


def is_reverse_palindrome(seq):
    return seq == revcomp(seq)

def palindrome_seq(dna_sequence,min_length=4,max_length=12):
    palindromes = []
    for length in range(min_length, max_length + 1):
        for i in range(len(dna_sequence) - length + 1):
            substring = dna_sequence[i:i+length]
            if is_reverse_palindrome(substring):
                palindromes.append((i + 1, length))
    return palindromes

def main():

    dna_sequence ='GCGCAGCCGTGAACAAGACTGAGTGCTTAGTTTGCTAATTGAGCTAGCGGAAACCTACTCTTTACTCACGTCAATTTATTTTCAGTGTACTGTGGACTAAATATCAAGTTCCATGCATTTCGGGTTGGCGCCCCATGATAGGACATCATTCATATGAATGCCGGCCCATCTTCATAAGCTGGTCCGGGGACATGGGCTTTTCACGCCGGTGCAACACCTGCCTGTTGAGTTAACGGAGCGCCGGGGGTAACGGACAACAGGATGGTGCCTGTTCGGGAGGCCCGACTCGGGCACAATTGATATTGTATGATCCGCTGTTATACCACACATGCGTAACAGTATGCTTTGCGAGGACCACCATGTTGGAGCGCAGTGTATGCGATTGTCTATCTACGGATAAGGAGCAACTCCCGACGGGCGTAAAGCCCACTAACCTTGAGACAACATCCGATCCTAGTGGACTTGTAGTACTACATTGGACCCGCATCTGTTCGAAATTATACCATAATCCGCTGACTAACTGCGCCAGAGAAAATGGTGCAAGCAACCTTAGTCGAGCGAGGGCACGCCGGAGCCTCATGCACCAGAGGCAGTAGAGCTCATTTACCCGCTCACCTCACGCCCGCACTCCTATGTTGTAGATTTTGACGAGCTACCTATTCCATTCGACTAGCTTTCCGCCTTCTCGCCGTCGGGTGCTCTGCGTGGCCAGCACTAGAGTTAAATGGCGCACGCAAACATAACAGTCGGCAGTTACCGCGGTTACCGGTAACTCACGTCCTTAGAATGCCCCAACTTGGTATATTCCCGGATTCACCCGCCGTCTGGCCGCACGACTTATTTGGCGGTGGTCAGCAGAAGTTTTGGCCAAGAACCTGGGGGATGCAATCAGTGACCACCTGCACAGCCTAACCTATTGCACTGGATACAAACAATAT'
    reverse_palindromes = palindrome_seq(dna_sequence)

    for start, length in reverse_palindromes:
        print(f" {start} {length}")

if __name__ == "__main__":
    main()

