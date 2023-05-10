from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np

#Fetch Mushroom Data Set
mushroom = fetch_openml(data_id=24)
mushroom_df = mushroom.data
mushroom_target = mushroom.target

mushroom_onehot_df = pd.concat([pd.get_dummies(mushroom_df[col], drop_first = False) for col in mushroom_df], axis=1, keys=mushroom_df.columns)
mushroom_onehot_target = pd.get_dummies(mushroom_target, drop_first=True)

mushroom_onehot_df.insert(loc=0, column='target', value=mushroom_onehot_target)
mushroom_onehot = mushroom_onehot_df.to_numpy(dtype=int)

for i in range(1, 4):
    mushroom_train, mushroom_test = train_test_split(mushroom_onehot, train_size=0.7, random_state=i)

    np.savetxt("mushroom_train_" + str(i) + ".csv", mushroom_train, delimiter=",", fmt='%i')
    np.savetxt("mushroom_test_" + str(i) + ".csv", mushroom_test, delimiter=",", fmt='%i')
