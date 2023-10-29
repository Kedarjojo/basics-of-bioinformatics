def revcomp(base):
    comp_bases = {'A':'T', 'G':'C','C':'G','T':'A'}
    rev_comp = ''
    
    reversed_string=inp_string[::-1]
    for base in reversed_string:
        
        rev_comp+=comp_bases[base]
    return rev_comp
   
inp_string=""
output=revcomp(inp_string)
print(output)
