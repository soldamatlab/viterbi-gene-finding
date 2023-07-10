import numpy as np

from lib.model.model import ALPHABET

NO_POINTER = -1


def viterbi(model, sequence):
    n_states = len(model.states)
    n_sequence = len(sequence)

    pointers = NO_POINTER * np.ones((n_sequence + 2, n_states), dtype=int)
    probabilities = -float('inf') * np.ones((n_sequence + 2, n_states))
    probabilities[0, model.index.start] = 0
    
    for c in range(n_sequence):
        for s in range(n_states):
            if not model.states[s].emission == sequence[c]:
                continue # it is initialised to -inf

            best_ancestor = NO_POINTER
            best_probability = -float('inf')
            for a in range(n_states):
                new_probability = probabilities[c, a] + model.states[a].probabilities[s]
                if new_probability > best_probability:
                    best_ancestor = a
                    best_probability = new_probability
            pointers[c+1, s] = best_ancestor
            probabilities[c+1, s] = best_probability

    best_ancestor = NO_POINTER
    best_probability = -float('inf')
    for a in range(n_states):
        new_probability = probabilities[n_sequence, a] + model.states[a].probabilities[model.index.end]
        if new_probability > best_probability:
            best_ancestor = a
            best_probability = new_probability
    pointers[n_sequence+1, model.index.end] = best_ancestor
    probabilities[n_sequence+1, model.index.end] = best_probability

    best_path = []
    state = np.argmax(probabilities[-1, :]) # will always be model.index.end
    for c in range(n_sequence+1, 0, -1):
        state = pointers[c, state]
        best_path.append(model.states[state].interpretation)
    best_path = best_path[0:-1]
    best_path.reverse()
    return best_path
    
