import numpy as np


def count_transitions(model, sequence, train_subsequences):
    for subsequence in train_subsequences:
        gene_border_transitions = []
        for gene in subsequence.genes:
            gene_border_transitions.append(gene[0])
            gene_border_transitions.append(gene[1] - 2)

        start = subsequence.start
        state = model.index.start
        for transition in gene_border_transitions:
            for i in range(start, transition):
                next_state = model.states[state].transitions.default[sequence[i]]
                model.states[state].probabilities[next_state] += 1
                state = next_state

            next_state = model.states[state].transitions.gene_border[sequence[transition]]
            model.states[state].probabilities[next_state] += 1
            state = next_state
            start = transition + 1

        for i in range(gene_border_transitions[-1] + 1, subsequence.end - 1):
            next_state = model.states[state].transitions.default[sequence[i]]
            model.states[state].probabilities[next_state] += 1
            state = next_state

        model.states[state].probabilities[model.index.end] += 1
    return model


def convert_counts_to_logprobs(model):
    for s in range(len(model.states)):
        if sum(model.states[s].probabilities) == 0:
            model.states[s].probabilities -= float('inf')
            continue
        model.states[s].probabilities /= sum(model.states[s].probabilities)
        np.seterr(divide='ignore')  # log(0) = -inf is intended here
        model.states[s].probabilities = np.log(model.states[s].probabilities)
        np.seterr(divide='warn')
    return model