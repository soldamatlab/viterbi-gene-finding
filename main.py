#!/usr/bin/python

from lib.model.model import Model
from lib.model.train import *
from lib.input import *
from lib.viterbi import viterbi
from lib.output import *

SEQUENCE_PATH = "sequenceOfEcoliStrainM54.txt"
TRAIN_PATH = "train.seqs"
TEST_PATH = "test.seqs"

if __name__ == '__main__':
    model = Model()

    sequence = read_sequence(SEQUENCE_PATH)
    train_subsequences = read_subsequences(TRAIN_PATH)

    # ___ train model ___________________________________________________________
    model = count_transitions(model, sequence, train_subsequences)
    model = convert_counts_to_logprobs(model)
    # ___________________________________________________________________________

    test_subsequences = read_subsequences(TEST_PATH)

    # ___ test model ____________________________________________________________
    predicted_subsequences = []
    for subsequence in test_subsequences:
        path = viterbi(model, sequence[subsequence.start:subsequence.end+1])
        predicted_subsequences.append(Subsequence(
            subsequence.start, subsequence.end, get_genes(path, subsequence.start)))
    # ___________________________________________________________________________

    with open("predictions.txt", "w") as file:
        for subsequence in predicted_subsequences:
            file.write("%d %d" % (subsequence.start+1, subsequence.end+1))
            for gene in subsequence.genes:
                file.write(" [%d, %d]" % (gene[0]+1, gene[1]+1))
            file.write("\n")

    n_predicted, n_overlapping, total_length, overlap_length, jaccard = get_accuracy(
        predicted_subsequences, test_subsequences)

    with open("accuracy.txt", "w") as file:
        file.write("%d\n" % (n_predicted))
        file.write("%d\n" % (n_overlapping))
        file.write("%d\n" % (total_length))
        file.write("%d\n" % (overlap_length))
        file.write("%f\n" % (jaccard))
