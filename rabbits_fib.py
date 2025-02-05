def rabbit_pairs(n, k):

    fib_sequence = [1, 1]
    
    for i in range(2, n):
        next_term = fib_sequence[i-1] + k * fib_sequence[i-2]
        fib_sequence.append(next_term)

    # Return the last term in the sequence
    return fib_sequence[-1]

print(rabbit_pairs(n,k))
