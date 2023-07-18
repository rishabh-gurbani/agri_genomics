import numpy as np
import pandas as pd

def nw(x, y, match=1, mismatch=1, gap=1):
    nx = len(x)
    ny = len(y)
    original_ref_length = ny

    # Truncate query sequence if it is longer than the reference sequence
    if nx > ny:
        x = x[:ny]
        nx = ny

    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, -nx * gap, nx + 1)
    F[0, :] = np.linspace(0, -ny * gap, ny + 1)

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1))
    P[:, 0] = 3
    P[0, :] = 4

    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if x[i] == y[j]:
                t[0] = F[i, j] + match
            else:
                t[0] = F[i, j] - mismatch
            t[1] = F[i, j + 1] - gap
            t[2] = F[i + 1, j] - gap
            tmax = np.max(t)
            F[i + 1, j + 1] = tmax
            if t[0] == tmax:
                P[i + 1, j + 1] += 2
            if t[1] == tmax:
                P[i + 1, j + 1] += 3
            if t[2] == tmax:
                P[i + 1, j + 1] += 4

    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:
            rx.append(x[i - 1])
            ry.append(y[j - 1])
            i -= 1
            j -= 1
        elif P[i, j] in [3, 5, 7, 9]:
            rx.append(x[i - 1])
            ry.append('-')
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j - 1])
            j -= 1

    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]

    # Trim the aligned sequences if they exceed the length of the original reference
    if len(rx) > original_ref_length:
        rx = rx[:original_ref_length]
        ry = ry[:original_ref_length]

    return rx, ry, F[nx, ny]


def normalise(query, reference):
    df = pd.read_csv("genotype_header.csv")
    df = df.iloc[:, 3:]
    if len(query) == len(reference):
        return query
    aligned, reference, score = nw(query, reference)
    for i in range(len(aligned)):
        if aligned[i] == "-":
            aligned = replace_character(aligned, i, df.iloc[5, i])
    return aligned


def replace_character(string, index, new_character):
    if index < 0 or index >= len(string):
        raise IndexError("Index is out of range")
    replaced_string = string[:index] + new_character + string[index + 1:]
    return replaced_string
