# -*- coding: utf-8 -*-
"""
Created on Sat May 15 18:59:48 2021

@author: x.zhou
"""

import tkinter as tk
from tkinter import filedialog
import tensorflow as tf
import numpy as np
import h5py
import os

root= tk.Tk()

canvas1 = tk.Canvas(root, width = 800, height = 240)
canvas1.pack()

def cnn():  
    data_path = filedialog.askopenfilename(filetypes=[("hdf5 file","*.hdf5")])
    file_path = os.path.split(data_path)[0]
    
    model_path = filedialog.askopenfilename(filetypes=[("h5 file","*.h5")])
    
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
    
    labelOut = open(file_path + "/APTdata_label.txt", "w+")
    for i in range(len(data_images)):
        labelOut.write("%i\n" % (data_value_predictions[i]))
    labelOut.close()
    
    label1 = tk.Label(root, text= 'Done! Please close the window!', fg='green', font=('helvetica', 12, 'bold'))
    canvas1.create_window(400, 180, window=label1)
    
button1 = tk.Button(text='Click Me',command=cnn, bg='brown',fg='white')
canvas1.create_window(400, 120, window=button1)
label1 = tk.Label(root, text= 'First choose .hdf5 file to import your experimental data, then choose the cnn model (.h5 file)', fg='green', font=('helvetica', 12, 'bold'))
canvas1.create_window(400, 60, window=label1)

root.mainloop()