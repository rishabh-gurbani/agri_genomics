import alignmentScore as aligner
def height_range(df, sequences,query_sequence, nk=5):
    import json

    # Read the JSON data from a file
    with open("json/clusters_dict.json", "r") as f:
        clusters_dict = json.loads(f.read())

    with open("json/sequence_cluster_dict.json", "r") as f:
        sequence_cluster_dict = json.loads(f.read())

    with open("json/sequences_dict.json", "r") as f:
        sequences_dict = json.loads(f.read())

    with open("json/unique_sequences_per_cluster.json", "r") as f:
        unique_sequences_per_cluster = json.loads(f.read())

    alignmentScore = dict()
    for k,v in sequences.items():
        alignmentScore[k] = aligner.alignSequencesGlobal(query_sequence, v)
    sorted_dict_by_value_desc = dict(sorted(alignmentScore.items(), key=lambda item: item[1], reverse=True))
    topk = list(sorted_dict_by_value_desc.keys())[:nk]
    clusters = []
    subpops = []
    for cultivar in topk:
        row = df.loc[df['Cultivar ID'] == cultivar]
        clusters.append(sequence_cluster_dict[cultivar])
        subpops.append(row['Subpopulation'].item())

    cluster_belong = sorted(clusters, key= lambda x : clusters.count(x), reverse=True)[0]
    subpop_belong = sorted(subpops, key = lambda x : subpops.count(x), reverse=True)[0]
    df_new = df[df['Cultivar ID'].isin(clusters_dict[str(cluster_belong)])]
    df_new = df_new[df_new['Subpopulation'] == subpop_belong]
    heights = df_new['Plant Height (cm)']
    # print("Height range : ", heights.min(), " - ", heights.min(),)
    return heights.min(),heights.max()