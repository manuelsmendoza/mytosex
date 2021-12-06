import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from keras.models import Sequential
from keras.layers import Dense

reference_data = pd.read_csv("/Users/manuelmendoza/Desktop/reference_data.tsv", delimiter="\t")
val = reference_data.iloc[:, 2:8]
sex = list(reference_data.loc[:, "sex"])
sex = [s.replace("Male", "1") for s in sex]
sex = [s.replace("Female", "0") for s in sex]
sex = [int(s) for s in sex]
sex = np.ravel(sex)

val_train, val_test, sex_train, sex_test = train_test_split(val, sex, test_size=1, random_state=40)

scaler = StandardScaler().fit(val_train)
val_train = scaler.transform(val_train)
val_test = scaler.transform(val_test)

model = Sequential()
model.add(Dense(9, activation='relu', input_dim=6))
model.add(Dense(6, activation='relu'))
model.add(Dense(3, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
model.fit(val_train, sex_train, epochs=20, batch_size=1, verbose=1)

model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
model.fit(val_train, sex_train, epochs=20, batch_size=1, verbose=1)

y_pred = model.predict(val_test)

