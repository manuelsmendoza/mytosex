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
from tensorflow import keras


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
#
# model = Sequential()
# model.add(Dense(6, activation='relu', input_dim=6))
# model.add(Dense(6, activation='relu'))
# model.add(Dense(3, activation='relu'))
# model.add(Dense(1, activation='sigmoid'))
# model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
#
# metrics_values = data_reference.loc[:, ["mtfcov", "mtmcov", "mtfmd", "mtmmd", "mtfgi", "mtmgi"]]
# samples_sex = data_reference.loc[:, "sex"]
# model.fit(metrics_values, samples_sex, epochs=250, batch_size=10, verbose=0)
# model.save(os.path.join(data_dir, "nn_model.h5"), save_format="h5")

print(tnow() + " INFO: Inferring the sex of the samples", file=sys.stdout)
model = load_model(os.path.join(model_dir, "nn_model.h5"))
align_metrics = pd.read_csv(
    os.path.join(data_dir, "align_stats.tsv"),
    delimiter="\t"
)

samples_name = list(align_metrics.loc[:, "sample"])
samples_info = align_metrics.loc[:, ["mtfcov", "mtmcov", "mtfmd", "mtmmd", "mtfgi", "mtmgi"]]
sex_prediction = model.predict(samples_info)
print(sex_prediction)
print(type(sex_prediction))
sex_prediction = [sex_round(x) for x in sex_prediction]
#sex_prediction = [int(x.round()) for x in sex_prediction]

results = {"sample": samples_name, "sex": sex_prediction}
results = pd.DataFrame.from_dict(results)
results["sex"].replace({0: "Female", 1: "Male"}, inplace=True)
results.to_csv(
    os.path.join(settings["output_dir"], "results.tsv"),
    sep="\t",
    index=False
)

plot_data = {
    "x_value": np.array(samples_info.loc[:, "mtmcov"] / samples_info.loc[:, "mtfcov"]),
    "y_value": np.array((samples_info.loc[:, "mtmgi"] * samples_info.loc[:, "mtmmd"]) /
                        (samples_info.loc[:, "mtfgi"] * samples_info.loc[:, "mtfmd"])),
    "sexpred": sex_prediction
}
plot_data = pd.DataFrame.from_dict(plot_data)
plot_data["sexpred"].replace({"Female": 0, "Male": 1}, inplace=True)

sns.scatterplot(
    data=plot_data,
    x="x_value",
    y="y_value",
    hue="sexpred")
# male_res = samples_info[samples_info["sex"] == "Male"]
# male_x = np.array(male_res.loc[:, "mtmcov"] / male_res.loc[:, "mtfcov"])
# male_y = np.array((male_res.loc[:, "mtmgi"] * male_res.loc[:, "mtmmd"]) /
#                   (male_res.loc[:, "mtfgi"] * male_res.loc[: , "mtfmd"]))
# male_c = np.repeat(1, male_res.shape[0])
# female_res = samples_info[samples_info["sex"] == "Female"]
#
#
#
# fig, ax = plt.subplots()
# ax.scatter(plot_data.loc[:, "x_value"], plot_data.loc[:, "y_value"], c=plot_data.loc[:, "sexpred"])
# yellow_patch = mpatches.Patch(color="yellow", label="Male")
# purple_patch = mpatches.Patch(color="Purple", label="Female")
# ax.legend(handles=[yellow_patch, purple_patch])
# plt.xlabel("Coverage")
# plt.ylabel("Sequencing depth")
plt.savefig(os.path.join(figs_dir, "clust_prediction.svg"))

open(os.path.join(settings["output_dir"], ".prediction.ok"), "w").close()
