from lib.model.model import NON_GENE


def get_genes(path, start):
    gene_borders = []
    previous = NON_GENE
    for i in range(len(path)):
        if path[i] != previous:
            if previous == NON_GENE:
                gene_borders.append(i)
            else:
                gene_borders.append(i-1)
        previous = path[i]

    genes = []
    for g in range(len(gene_borders)):
        if not g % 2:
            continue
        genes.append([gene_borders[g-1]+start, gene_borders[g]+start])

    return genes


def get_accuracy(predicted_subsequences, test_subsequences):
    n_predicted = 0
    n_overlapping = 0
    total_length = 0
    overlap_length = 0

    for s in range(len(predicted_subsequences)):
        n_predicted += len(predicted_subsequences[s].genes)

        for gene in predicted_subsequences[s].genes:
            overlap = False
            for true_gene in test_subsequences[s].genes:
                ol = min([gene[1], true_gene[1]]) - max([gene[0], true_gene[0]]) + 1
                ol = max(ol, 0)
                overlap_length += ol

                if ol > 0:
                    overlap = True
            if overlap:
                n_overlapping += 1
                
        for gene in predicted_subsequences[s].genes:
            total_length += gene[1] - gene[0] + 1

    test_total_length = 0
    for subsequence in test_subsequences:
        for gene in subsequence.genes:
            test_total_length += gene[1] - gene[0] + 1

    jaccard = (overlap_length) / (total_length + test_total_length - overlap_length)

    return n_predicted, n_overlapping, total_length, overlap_length, jaccard
