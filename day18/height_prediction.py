import pickle
from subPopulation import subPopulation
from sequence_to_ohe import sequence_to_ohe


def predict(sequence):

    with open('model.pkl', 'rb') as f:
        model = pickle.load(f)

    normalised_seq, subPop = subPopulation(sequence)
    ohe = sequence_to_ohe(normalised_seq, subPop)
    prediction = model.predict(ohe)

    return prediction

