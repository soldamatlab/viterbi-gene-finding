import numpy as np

from lib.model.index import Index

ALPHABET = ['A', 'C', 'G', 'T']
REGULARIZATION_CONST = 1

GENE = 1
NON_GENE = 0


class Interpretation:
    START = NON_GENE
    END = NON_GENE
    NONCODING = NON_GENE
    START_CODON = GENE
    INTERNAL_CODONS = GENE
    STOP_CODON = GENE


class Transitions:
    def __init__(self):
        self.default = dict.fromkeys(ALPHABET, None)
        self.gene_border = dict()


class State:
    def __init__(self, interpretation, emission):
        self.interpretation = interpretation
        self.emission = emission
        self.probabilities = None
        self.transitions = Transitions()

    def init_transitions(self, n_states):
        self.probabilities = np.zeros(n_states)

    def add_transition(self, to, default=None, gene_border=None):
        self.probabilities[to] = REGULARIZATION_CONST
        if default is not None:
            self.transitions.default[default] = to
        if gene_border is not None:
            self.transitions.gene_border[gene_border] = to


class Model:
    def __init__(self):
        self.states = []
        self.index = Index()

        self.init_states()
        self.init_transitions()

    def add_state(self, index_where, index_attribute, interpretation, emission):
        setattr(index_where, index_attribute, len(self.states))
        self.states.append(State(interpretation, emission))

    def init_states(self):
        # start
        self.add_state(self.index, 'start', Interpretation.START, '')

        # noncoding
        for char in ALPHABET:
            self.add_state(self.index.noncoding, char,
                           Interpretation.NONCODING, char)

        # start codon
        for char in ALPHABET:
            self.add_state(self.index.start_codon.first, char,
                           Interpretation.START_CODON, char)
        self.add_state(self.index.start_codon, 'second',
                       Interpretation.START_CODON, 'T')
        self.add_state(self.index.start_codon, 'third',
                       Interpretation.START_CODON, 'G')

        # internal codons
        for char in ALPHABET:
            self.add_state(self.index.internal_codons.first,
                           char, Interpretation.INTERNAL_CODONS, char)
        for char in ALPHABET:
            self.add_state(self.index.internal_codons.second,
                           char, Interpretation.INTERNAL_CODONS, char)
        for char in ALPHABET:
            self.add_state(self.index.internal_codons.third,
                           char, Interpretation.INTERNAL_CODONS, char)

        # stop codon
        self.add_state(self.index.stop_codon, 'first',
                       Interpretation.STOP_CODON, 'T')
        for char in ['A', 'G']:
            self.add_state(self.index.stop_codon.second, char,
                           Interpretation.STOP_CODON, char)
        for char in ['A', 'G']:
            self.add_state(self.index.stop_codon.third, char,
                           Interpretation.STOP_CODON, char)

        # end
        self.add_state(self.index, 'end', Interpretation.END, '')

    def init_transitions(self):
        for s in range(len(self.states)):
            self.states[s].init_transitions(len(self.states))
        self.init_internal_transitions()
        self.init_external_transitions()

    def init_internal_transitions(self):
        # noncoding
        for s in ALPHABET:
            for t in ALPHABET:
                self.states[getattr(self.index.noncoding, s)].add_transition(
                    getattr(self.index.noncoding, t), default=t)

        # start codon
        for char in ALPHABET:
            self.states[getattr(self.index.start_codon.first, char)].add_transition(
                getattr(self.index.start_codon, 'second'), default='T')
        self.states[getattr(self.index.start_codon, 'second')].add_transition(
            getattr(self.index.start_codon, 'third'), default='G')

        # internal codons
        for s in ALPHABET:
            for t in ALPHABET:
                self.states[getattr(self.index.internal_codons.first, s)].add_transition(
                    getattr(self.index.internal_codons.second, t), default=t)
                self.states[getattr(self.index.internal_codons.second, s)].add_transition(
                    getattr(self.index.internal_codons.third, t), default=t)

        # stop codon
        for char in ['A', 'G']:
            self.states[self.index.stop_codon.first].add_transition(
                getattr(self.index.stop_codon.second, char), default=char)
        for s in ['A', 'G']:
            for t in ['A', 'G']:
                self.states[getattr(self.index.stop_codon.second, s)].add_transition(
                    getattr(self.index.stop_codon.third, t), default=t)

    def init_external_transitions(self):
        # start -> noncoding
        for char in ALPHABET:
            self.states[self.index.start].add_transition(
                getattr(self.index.noncoding, char), default=char)

        # noncoding -> start codon
        for noncoding in ALPHABET:
            for start_codon in ALPHABET:
                self.states[getattr(self.index.noncoding, noncoding)].add_transition(
                    getattr(self.index.start_codon.first, start_codon), gene_border=start_codon)

        # start codon -> internal codons
        for char in ALPHABET:
            self.states[self.index.start_codon.third].add_transition(
                getattr(self.index.internal_codons.first, char), default=char)

        # internal codons -> internal codons
        for third in ALPHABET:
            for first in ALPHABET:
                self.states[getattr(self.index.internal_codons.third, third)].add_transition(
                    getattr(self.index.internal_codons.first, first), default=first)

        # internal codons -> stop codon
        for char in ALPHABET:
            self.states[getattr(self.index.internal_codons.third, char)].add_transition(
                self.index.stop_codon.first, gene_border='T')

        # stop codon -> noncoding
        for stop_codon in ['A', 'G']:
            for noncoding in ALPHABET:
                self.states[getattr(self.index.stop_codon.third, stop_codon)].add_transition(
                    getattr(self.index.noncoding, noncoding), default=noncoding)

        # noncoding -> end
        for char in ALPHABET:
            self.states[getattr(self.index.noncoding, char)].add_transition(
                self.index.end)
