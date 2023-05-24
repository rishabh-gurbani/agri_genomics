from enum import IntEnum
import numpy as np
def printMatrix(mat):
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            print(mat[i][j], end=" ")
        print()

def constructMatrix(str1, str2):
    list = [[0 for i in range(len(str1)+1)] for j in range(len(str2)+1)]
    return list

def printStringMatrix(matrix, str1, str2):
    print("", end = "  ")
    print("_", end = " ")
    for i in range(len(str1)):
        print(str1[i], end = " ")
    print()
    for i in range(len(matrix)):
        if i == 0 :
            print("_", end = " ")
        else:
            print(str2[i-1], end = " ")
        for j in range(len(matrix[i])):
            print(matrix[i][j], end = " ")
        print()

def initMatrix(matrix, gapPenalty):
    x=0
    for i in range(len(matrix)):
        matrix[i][0] = x * gapPenalty
        x=x+1
    x=0
    for i in range(len(matrix[0])):
        matrix[0][i] = x * gapPenalty
        x=x+1
    return matrix

def getMinPenalty(matrix, i, j, str1, str2,matchPenalty, mismatchPenalty, gapPenalty ):
    if str2[i-1] == str1[j-1]:
        currPenalty = matchPenalty
    else:
        currPenalty = mismatchPenalty

    top = matrix[i-1][j] + gapPenalty
    left = matrix[i][j-1] + gapPenalty
    diag = matrix[i-1][j-1] + currPenalty
    maxE = max(top, left, diag)
    currPenalty = maxE

    return currPenalty

def fillMatrix(matrix, str1, str2, matchPenalty, mismatchPenalty, gapPenalty):
    for i in range(1, len(matrix)):
        for j in range(1, len(matrix[0])):
            matrix[i][j] = getMinPenalty(matrix, i, j, str1, str2, matchPenalty, mismatchPenalty, gapPenalty)


def getPrevious(matrix, i, j,matchPenalty, mismatchPenalty, gapPenalty ):

    score = matrix[i][j]
    top = matrix[i-1][j]
    left = matrix[i][j-1]
    diag = matrix[i-1][j-1]

    maxPrevIndexI = 0
    maxPrevIndexJ = 0
    maxPrevDir = 0
    if diag + matchPenalty == score or diag + mismatchPenalty == score:
        maxPrevIndexI = i-1
        maxPrevIndexJ = j-1
        maxPrevDir = 0
    elif top + gapPenalty == score:
        maxPrevIndexI = i-1
        maxPrevIndexJ = j
        maxPrevDir = 1
    elif left + gapPenalty == score:
        maxPrevIndexI = i
        maxPrevIndexJ = j-1
        maxPrevDir = 2

    return maxPrevDir, maxPrevIndexI, maxPrevIndexJ

def score(matrix):
    return matrix[-1][-1]

def backTrack(matrix, matchPenalty, mismatchPenalty, gapPenalty):

    score = 0
    directions = []

    i = len(matrix)-1
    j = len(matrix[0])-1

    while i>=0 and j>=0:
        score += matrix[i][j]
        if i == 0 or j == 0:
            break
        prev = getPrevious(matrix, i, j, matchPenalty, mismatchPenalty, gapPenalty)
        directions.insert(0, prev[0])
        i = prev[1]
        j = prev[2]
    return directions, score


def alignSequencesGlobal(sequence1, sequence2, matchPenalty = 1, mismatchPenalty = -1, gapPenalty = -2):
    mat = constructMatrix(sequence1, sequence2)
    matrix = initMatrix(mat, gapPenalty)
    fillMatrix(matrix, sequence1, sequence2, matchPenalty, mismatchPenalty, gapPenalty)
    return matrix[-1][-1]
    # dir = backTrack(matrix, matchPenalty, mismatchPenalty, gapPenalty)
    # return dir[1]


# Assigning the constants for the scores
class Score(IntEnum):
    MATCH = 1
    MISMATCH = -1
    GAP = -1

# Assigning the constant values for the traceback
class Trace(IntEnum):
    STOP = 0
    LEFT = 1
    UP = 2
    DIAGONAL = 3

def alignSequencesLocal(seq1, seq2):
    # Generating the empty matrices for storing scores and tracing
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)
    tracing_matrix = np.zeros(shape=(row, col), dtype=int)

    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)

    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = Score.MATCH if seq1[i - 1] == seq2[j - 1] else Score.MISMATCH
            diagonal_score = matrix[i - 1, j - 1] + match_value

            # Calculating the vertical gap score
            vertical_score = matrix[i - 1, j] + Score.GAP

            # Calculating the horizontal gap score
            horizontal_score = matrix[i, j - 1] + Score.GAP

            # Taking the highest score
            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)

            # Tracking where the cell's value is coming from
            if matrix[i, j] == 0:
                tracing_matrix[i, j] = Trace.STOP

            elif matrix[i, j] == horizontal_score:
                tracing_matrix[i, j] = Trace.LEFT

            elif matrix[i, j] == vertical_score:
                tracing_matrix[i, j] = Trace.UP

            elif matrix[i, j] == diagonal_score:
                tracing_matrix[i, j] = Trace.DIAGONAL

            # Tracking the cell with the maximum score
            if matrix[i, j] >= max_score:
                max_index = (i,j)
                max_score = matrix[i, j]

    # Initialising the variables for tracing
    aligned_seq1 = ""
    aligned_seq2 = ""
    current_aligned_seq1 = ""
    current_aligned_seq2 = ""
    (max_i, max_j) = max_index

    # Tracing and computing the pathway with the local alignment
    while tracing_matrix[max_i, max_j] != Trace.STOP:
        if tracing_matrix[max_i, max_j] == Trace.DIAGONAL:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1

        elif tracing_matrix[max_i, max_j] == Trace.UP:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'
            max_i = max_i - 1

        elif tracing_matrix[max_i, max_j] == Trace.LEFT:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]
            max_j = max_j - 1

        aligned_seq1 = aligned_seq1 + current_aligned_seq1
        aligned_seq2 = aligned_seq2 + current_aligned_seq2

    # Reversing the order of the sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return aligned_seq1, aligned_seq2, max_score
