import pandas as pd
import pickle
import os

MODEL_PATH = os.path.join(os.getcwd(), 'micro_visualize', 'models', 'rf_resistance.pkl')

def load_model():
    with open(MODEL_PATH, 'rb') as f:
        model = pickle.load(f)
    return model

def predict_resistance(df):
    model = load_model()
    features = df.select_dtypes(include='number')
    predictions = model.predict(features)
    df['Predicted_Resistance'] = predictions
    return df
