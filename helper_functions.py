from Bio.Align import substitution_matrices

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    # Initialise lengths of sequences and matrices
    n, m = len(seq1), len(seq2)
    scores = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    arrows = [[None for _ in range(m + 1)] for _ in range(n + 1)]
  
    for i in range(1, n + 1):
        scores[i][0] = scores[i-1][0] + scoring_function(seq1[i - 1], "-")
        arrows[i][0] = "up"

    for j in range(1, m + 1):
        scores[0][j] = scores[0][j-1] + scoring_function("-", seq2[j - 1])
        arrows[0][j] = "left"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Fill in scores matrix with scores using Needleman-Wunsch algorithm
            match = scores[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = scores[i - 1][j] + scoring_function(seq1[i - 1], "-")
            left = scores[i][j - 1] + scoring_function("-", seq2[j - 1])
            scores[i][j] = max(match, up, left)

            # Fill in arrows matrix with directions using results from algorithm
            if scores[i][j] == match:
                arrows[i][j] = "diagonal"
            elif scores[i][j] == up:
                arrows[i][j] = "up"
            else:
                arrows[i][j] = "left"

    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n, m

    # Find full path back to start of matrix to get aligned sequences
    while i > 0 or j > 0:
        if arrows[i][j] == "diagonal":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif arrows[i][j] == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse directions to get aligned sequences
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    aligned_seq1_combined = "".join(aligned_seq1)
    aligned_seq2_combined = "".join(aligned_seq2)

    return aligned_seq1_combined, aligned_seq2_combined, float(scores[n][m])

def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    # Initialise lengths of sequences and matrices
    n, m = len(seq1), len(seq2)
    scores = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    arrows = [[None for _ in range(m + 1)] for _ in range(n + 1)]
  
    for i in range(1, n + 1):
        scores[i][0] = scores[i-1][0] + scoring_function(seq1[i - 1], "-")
        arrows[i][0] = "up"

    for j in range(1, m + 1):
        scores[0][j] = scores[0][j-1] + scoring_function("-", seq2[j - 1])
        arrows[0][j] = "left"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Fill in scores matrix with scores using Needleman-Wunsch algorithm
            match = scores[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = scores[i - 1][j] + scoring_function(seq1[i - 1], "-")
            left = scores[i][j - 1] + scoring_function("-", seq2[j - 1])
            scores[i][j] = max(match, up, left)

            # Fill in arrows matrix with directions using results from algorithm
            if scores[i][j] == match:
                arrows[i][j] = "diagonal"
            elif scores[i][j] == up:
                arrows[i][j] = "up"
            else:
                arrows[i][j] = "left"

    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n, m

    # Find full path back to start of matrix to get aligned sequences
    while i > 0 or j > 0:
        if arrows[i][j] == "diagonal":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif arrows[i][j] == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse directions to get aligned sequences
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    aligned_seq1_combined = "".join(aligned_seq1)
    aligned_seq2_combined = "".join(aligned_seq2)

    return aligned_seq1_combined, aligned_seq2_combined, float(scores[n][m])


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

# Scoring matrix function using BLOSUM62 and a gap penalty of -8
def scoring_function(aa_i,aa_j):
    blosum62 = substitution_matrices.load("BLOSUM62")
    gap_penalty = -8

    if aa_i == "-" or aa_j == "-":
        score = gap_penalty
    else:
        score = blosum62[aa_i, aa_j]

    return (score)
