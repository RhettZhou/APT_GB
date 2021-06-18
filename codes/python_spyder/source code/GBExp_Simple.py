# -*- coding: utf-8 -*-
"""

This file is a ealier version of cnn.py. 

The purpose of this file is to export labels for each of the individual images 
that were saved in the .hdf5 file. 
The trained CNN model was saved in the file 'model.h5'. 
This .h5 should should be in the same folder as this .py file. 

Created on April 17 18:59:48 2020

@author: x.zhou
"""
from tkinter import filedialog
import tensorflow as tf
import numpy as np
import h5py
import os

### Process the data

data_path = filedialog.askopenfilename(filetypes=[("hdf5 file","*.hdf5")])
file_path = os.path.split(data_path)[0]

model_path = filedialog.askopenfilename()

hf = h5py.File(data_path, 'r')
data_images = []
for k in hf.keys():
    data = hf[k][:]
    data_images.append(data)

data_images_new = np.expand_dims(data_images, axis=-1)
print(data_images_new.shape)

new_model = tf.keras.models.load_model(model_path)
new_model.summary()
new_model.load_weights(model_path)

data_predictions = new_model.predict(data_images_new)

data_value_predictions = np.zeros(len(data_images))

for i in range(len(data_images)):
    data_value_predictions[i] = np.argmax(data_predictions[i])

labelOut = open(file_path + "/APTdata_label.txt", "w+")
for i in range(len(data_images)):
    labelOut.write("%i\n" % (data_value_predictions[i]))
labelOut.close()
