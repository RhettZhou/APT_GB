# -*- coding: utf-8 -*-
"""

This file is to plot the confusion matrix for the testing data. 
One needs to select the folder that contains the testing images and lables.
One also needs to select a trained CNN model. 

Created on Thu Apr 16 20:35:30 2020

@author: x.zhou
"""
from tkinter import filedialog
import tensorflow as tf
import numpy as np

import os
import cv2
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

### Process the data
def load_image(file_path):
    return cv2.imread(file_path,0)

file_path = filedialog.askdirectory()
model_path = filedialog.askopenfilename()

test_path = file_path + "/test/"
test_image_files = os.listdir(test_path)
test_images_length = len(test_image_files)
test_images = [load_image(test_path + str(file+1) + '.tif') for file in range(test_images_length)]
test_labels_file = open(file_path + "/test.txt", "r")
test_labels = test_labels_file.readlines()
test_labels_file.close()
for i in range(len(test_images)):
   test_images[i] = test_images[i] / 255.0
test_images_new = np.expand_dims(test_images, axis=-1)
test_labels = np.array(test_labels)
print(test_images_new.shape, test_labels.shape)   

new_model = tf.keras.models.load_model(model_path)
new_model.summary()
new_model.load_weights(model_path)

test_evaluate = new_model.evaluate(test_images_new, test_labels)
test_predictions = new_model.predict(test_images_new)

### Visulazation  Test images
cols = 5
rows = np.ceil(len(test_images)/cols)
fig = plt.figure()
fig.set_size_inches(cols*2, rows*2)
value_label = np.zeros(len(test_images))
value_predictions = np.zeros(len(test_images))
ii = -1;
for i in range(len(test_images)):
    value_label[i] = int(test_labels[i])
    value_predictions[i] = np.argmax(test_predictions[i])
    check_predictions = "True" if value_label[i] == value_predictions[i] else "False"
    text_predictions = "InGrain" if value_predictions[i] == 0 else 'GB' if value_predictions[i] == 1 else 'Tri'
    if check_predictions == "False":
        ii = ii + 1
        plt.subplot(rows, cols, ii+1)
        plt.imshow(test_images[i])
        plt.title([text_predictions + '_' + check_predictions])
        plt.axis('off')

def plot_confusion_matrix(df_confusion, title='Confusion matrix', cmap=plt.cm.gray_r):
    plt.matshow(df_confusion, cmap=cmap) # imshow
    plt.colorbar()
    plt.ylabel('Predictions')
    plt.xlabel('Actual values')

test_confusion = confusion_matrix(value_label, value_predictions)
fig = plt.figure()
plot_confusion_matrix(test_confusion)
