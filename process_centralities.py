import pickle
import pandas as pd

if __name__ == '__main__':
    
    with open('aws-downloads/centralities.pickle', 'rb') as handle:
        data = pickle.load(handle)

        centralidades_perturbadas = data['centralidades_perturbadas']
        baseline = data['baseline']
        
        print(f"Read centralidades_perturbadas of length: {len(centralidades_perturbadas)}")
        print(f"Read baseline of length: {len(baseline)}")
        baseline.to_csv('1aaaa.csv')
        