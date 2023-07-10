UNDEFINED_INDEX = -1


class Quadruple:
    A = UNDEFINED_INDEX
    C = UNDEFINED_INDEX
    G = UNDEFINED_INDEX
    T = UNDEFINED_INDEX


class StopCodonConsecutive:
    A = UNDEFINED_INDEX
    G = UNDEFINED_INDEX


class StartCodon:
    first = Quadruple()
    second = UNDEFINED_INDEX
    third = UNDEFINED_INDEX


class InternalCodons:
    first = Quadruple()
    second = Quadruple()
    third = Quadruple()


class StopCodon:
    first = UNDEFINED_INDEX
    second = StopCodonConsecutive()
    third = StopCodonConsecutive()



class Index:
    """
    Stores indexes of states in Model.
    """
    start = UNDEFINED_INDEX
    noncoding = Quadruple()
    start_codon = StartCodon()
    internal_codons = InternalCodons()
    stop_codon = StopCodon()
    end = UNDEFINED_INDEX

