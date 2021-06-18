# -*- coding: utf-8 -*-
"""
The purpose of this file to compare the accuary of the trained CNN models.

Created on Thu Apr 16 19:20:21 2020

@author: x.zhou
"""

from tkinter import filedialog
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

folder_path = filedialog.askdirectory()
name_end = "_test-tag-accuracy.csv"

times = 5
maxData = 20

dense_layers = [3]    # 3layers is the best [1,2,3]
layer_sizes = [64] # 64 is better [32,64]
conv_layers = [1]  # 1 is better [1,2] very matter
batch_sizes = [32] # 32 is better [64,32] does not matter
lrs = [0.001]  # does not matter

record = np.zeros((maxData,times))
x = np.zeros(maxData)
row_num = 99

for no in range(1, maxData*times+1):
  noi = math.ceil(no/times)
  subi = no%times
  for batch_size in batch_sizes:
   for lr in lrs:
    for dense_layer in dense_layers:
     for layer_size in layer_sizes:
      for conv_layer in conv_layers:
        name = "-{}c-{}n-{}d-{}b-{}l".format(conv_layer,layer_size,dense_layer,batch_size,lr)
        name = str("{:0>5d}".format(noi*1000)) + "-" + str(subi) + name + name_end
        file_path = folder_path + "/run-S" + name
        mydataset = pd.read_csv(file_path)
        record[noi-1][subi] = mydataset.iloc[row_num][2]
        x[noi-1] = noi*1000
        
acc_mean = np.mean(record, axis=1)
acc_std = np.std(record, axis=1)
plt.plot(x, acc_mean, 'k-')
plt.fill_between(x, acc_mean-acc_std, acc_mean+acc_std)
plt.show()
#plt.errorbar(x, acc_mean, acc_std, linestyle='None', marker='^')
#plt.show()

