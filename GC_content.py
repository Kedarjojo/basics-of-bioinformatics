

def gc_content(sequence):
    gc_count = 0
    for base in sequence:
        if base == 'G' or base == 'C':
            gc_count += 1

    gc_content = (gc_count / len(sequence)) * 100
    return gc_content

fasta_file = "rosalind_gc.txt"
highest_gc_percent = 0
highest_gc_sequence = ""

with open(fasta_file, 'r') as file:
    sequence = ""
    for line in file:
        if line.startswith(">"):
            if sequence:
                current_gc_percent = gc_content(sequence)
                if current_gc_percent > highest_gc_percent:
                    highest_gc_percent = current_gc_percent
                    highest_gc_sequence = sequence
                sequence = ""
        else:
            sequence += line.strip()


print(f"Highest GC Content: {highest_gc_percent:.6f}%")
print("Sequence with Highest GC Content:")
print(highest_gc_sequence)






