# -*- coding: utf-8 -*-
"""
Created on Sat May 15 18:59:48 2021

@author: x.zhou
"""

import tensorflow as tf
import numpy as np
import h5py
import os
import argparse
 
parser = argparse.ArgumentParser(description='Program Arguments')
parser.add_argument('dp', help="Input filename, must be .hdf5")
parser.add_argument('--mp', default='cnn/model.h5')
opts = parser.parse_args()

data_path = opts.dp
file_path = os.path.split(data_path)[0]
file_name = os.path.split(data_path)[1]

model_path = opts.mp

hf = h5py.File(data_path, 'r')
data_images = []
for k in hf.keys():
    data = hf[k][:]
    data_images.append(data)

data_images_new = np.expand_dims(data_images, axis=-1)
#print(data_images_new.shape)

new_model = tf.keras.models.load_model(model_path)
new_model.summary()
new_model.load_weights(model_path)

data_predictions = new_model.predict(data_images_new)

data_value_predictions = np.zeros(len(data_images))

for i in range(len(data_images)):
    data_value_predictions[i] = np.argmax(data_predictions[i])

labelOut = open(file_path + "/" + file_name[:-5] + ".txt", "w+")
for i in range(len(data_images)):
    labelOut.write("%i\n" % (data_value_predictions[i]))
labelOut.close()