import numpy as np
import pandas as pd

def sequence_to_ohe(sequence, subpop):
    ohe_mapping = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}
    ohe_sequence = []

    for char in sequence:
        if char in ohe_mapping:
            ohe_sequence.append(ohe_mapping[char])
        else:
            ohe_sequence.append([0, 0, 0, 0])

    label_mapping = {
        'Aus': 0,
        'Indica': 1,
        'Intermediate': 2,
        'Japonica': 3,
        'VI/Aromatic': 4
    }

    encoded_label = label_mapping[subpop]

    ohe = np.append(np.array(ohe_sequence).flatten(), encoded_label)

    df = pd.read_csv("encoded.csv").iloc[:, :97]

    feature_names = df.columns.tolist()

    input_array = np.reshape(ohe, (-1, len(feature_names)))

    input_array_named = pd.DataFrame(input_array, columns=feature_names)

    return input_array_named




