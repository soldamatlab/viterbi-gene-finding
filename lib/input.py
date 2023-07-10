import re


class Subsequence:
    def __init__(self, start, end, genes=[]):
        self.start = start
        self.end = end
        self.genes = genes


def read_subsequences(path):
    subsequences = []
    with open(path, "r") as file:
        train = file.read().splitlines()
        for l in range(len(train)):
            subsequence = re.split(r'\D+', train[l])
            if subsequence[-1] == '':
                subsequence.pop()
            subsequence = list(map(int, subsequence))

            genes = []
            for i in range(3, len(subsequence)):
                if i % 2:
                    genes.append([subsequence[i-1]-1, subsequence[i]-1])
            subsequences.append(
                Subsequence(subsequence[0]-1, subsequence[1]-1, genes))
    return subsequences


def read_sequence(path):
    with open(path, "r") as file:
        sequence = file.read().replace('\n', '').upper()
    return sequence
