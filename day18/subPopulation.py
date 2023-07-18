import global_alignment as ga
from Bio import SeqIO
import pandas as pd
import json
from collections import Counter
import normalise_query


def subPopulation(query_sequence, nk = 5):
    sequences = {}
    for record in SeqIO.parse("sequences_24.fasta", "fasta"):
        sequences[record.id] = str(record.seq)

    df = pd.read_csv("mapping_24.csv")
    df['Subpopulation'] = df['Subpopulation'].astype("string")

    with open("json/clusters_dict.json", "r") as f:
        clusters_dict = json.loads(f.read())

    with open("json/sequence_cluster_dict.json", "r") as f:
        sequence_cluster_dict = json.loads(f.read())

    with open("json/sequences_dict.json", "r") as f:
        sequences_dict = json.loads(f.read())

    alignmentScore = dict()
    for k,v in sequences.items():
        query, actual, score = ga.nw(query_sequence, v)
        alignmentScore[k] = score

    sorted_dict_by_value_desc = dict(sorted(alignmentScore.items(),
                                            key=lambda item: item[1], reverse=True))
    topk = list(sorted_dict_by_value_desc.keys())[:nk]

    top_match = sequences_dict[topk[0]]
    normalised = normalise_query.normalise(query_sequence, top_match)

    clusters = []
    for cultivar in topk:
        row = df.loc[df['Cultivar ID'] == cultivar]
        clusters.append(sequence_cluster_dict[cultivar])

    cluster_belong = sorted(clusters, key= lambda x : clusters.count(x), reverse=True)[0]
    subpops = []
    for cultivar in clusters_dict[str(cluster_belong)]:
        x = df.loc[df['Cultivar ID']==cultivar]['Subpopulation'].values[0]
        subpops.append(x)

    # Count the occurrences of each element
    counter = Counter(subpops)

    # Find the value with the highest frequency
    most_common_value = counter.most_common(1)[0][0]
    return normalised, most_common_value
