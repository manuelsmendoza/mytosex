import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
import tensorflow as tf
from src.func_analysis import *
from src.func_setup import load_settings
from src.func_util import tnow, pass_file
from sklearn.preprocessing import StandardScaler
from keras.models import Sequential, save_model, load_model
from keras.layers import Dense
from keras.utils.vis_utils import plot_model
from tensorflow import keras
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import umap.umap_ as umap


settings = load_settings(os.getenv("MYTOSEX_SETTINGS"))
tmp_dir = os.path.join(settings["output_dir"], "tmp")
data_dir = os.path.join(settings["output_dir"], "data")
figs_dir = os.path.join(settings["output_dir"], "figures")
model_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "model")

# print(tnow() + " INFO: Building the neural network", file=sys.stdout)
# data_reference = pd.read_csv(
#     os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "model", "data.tsv"),
#     sep="\t"
# )
# data_reference["sex"].replace({"Female": 0, "Male": 1}, inplace=True)

# model = Sequential()
# model.add(Dense(6, activation='relu', input_dim=6))
# model.add(Dense(3, activation='relu'))
# model.add(Dense(1, activation='sigmoid'))
# model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
#
# data_reference = pd.read_csv(os.path.join(model_dir, "data.tsv"), delimiter="\t")
# data_reference["sex"].replace({"Female": 0, "Male": 1}, inplace=True)
# metrics_values = data_reference.loc[:, ["mtfcov", "mtmcov", "mtfmd", "mtmmd", "mtfgi", "mtmgi"]]
# samples_sex = data_reference.loc[:, "sex"]
# model.fit(metrics_values, samples_sex, epochs=250, batch_size=10)
# model.save(os.path.join(data_dir, "nn_model.h5"), save_format="h5")

print(tnow() + " INFO: Inferring the sex of the samples", file=sys.stdout)
model = load_model(os.path.join(model_dir, "med_model"))
align_metrics = pd.read_csv(
    os.path.join(data_dir, "align_stats.tsv"),
    delimiter="\t"
)

samples_name = list(align_metrics.loc[:, "sample"])
samples_info = align_metrics.loc[:, ["mtfcov", "mtmcov", "mtfmd", "mtmmd", "mtfgi", "mtmgi"]]
sex_prediction = model.predict(samples_info)
sex_prediction = np.array([x[0] for x in sex_prediction])
prediction_val = []
for pred in sex_prediction:
    if pred is np.NaN:
        prediction_val.append(pred)
    else:
        prediction_val.append(pred.round())

results = {"sample": samples_name, "sex": prediction_val}
results = pd.DataFrame.from_dict(results)
results["sex"].replace({0: "Female", 1: "Male"}, inplace=True)

# Add sexual index
results["sindex"] = align_metrics.loc[:, "mtmcov"] / align_metrics.loc[:, "mtfcov"]

# Dimension reduction
reducer = umap.UMAP()
dim_data = align_metrics[["mtfcov", "mtmcov", "mtfmd", "mtmmd", "mtfgi", "mtmgi"]].values
scaled_dim_data = StandardScaler().fit_transform(dim_data)
embedding = reducer.fit_transform(scaled_dim_data)
embedding = pd.DataFrame(embedding)
results["umapx"] = embedding.iloc[:, 0]
results["umapy"] = embedding.iloc[:, 1]

results.to_csv(
    os.path.join(settings["output_dir"], "results.tsv"),
    sep="\t",
    index=False
)

open(os.path.join(settings["output_dir"], ".prediction.ok"), "w").close()
