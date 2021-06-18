# -*- coding: utf-8 -*-
"""

After generating the images and labels from the Matlab script, one could this 
file to to train CNN models with different hyper parameters. 

Created on Sun Apr 12 00:16:35 2020

@author: x.zhou
"""

import tensorflow as tf
import math
import numpy as np
import os
import cv2
from tkinter import filedialog

#def makeWindowsCmdPath(path):
#    return '\"' + str(path) + '\"'
def load_image(file_path):
    return cv2.imread(file_path,0)

initial = 96
times = 5
maxData = 20
EPOCHS = 100

folder_path = filedialog.askdirectory()
split_type = "S"
folder_path = folder_path + "/" + split_type
model_path = folder_path + "/models"
try:
    os.stat(model_path)
except:
    os.mkdir(model_path) 

logs_base_dir =  folder_path + "/logs"
logs_base_dir = logs_base_dir.replace('/', '\\')
try:
    os.stat(logs_base_dir)
except:
    os.mkdir(logs_base_dir) 


dense_layers = [3]    # 3layers is the best [1,2,3]
dense_layer_sizes = [64,128,512]  # 
layer_sizes = [64] # 64 is better [32,64]
conv_layers = [1]  # 1 is better [1,2] very matter
batch_sizes = [32] # 32 is better [64,32] does not matter
lrs = [0.001]  # does not matter



for no in range(initial, maxData*times+1):
# if (math.ceil(no/times))%2 == 0:
    noi = math.ceil(no/times)
#    file_path = folder_path + "/" + split_type + str("{:0>5d}".format(noi*1000))
    file_path = folder_path + "/" + split_type + str(noi*1000)
    print(noi*1000,no%times)

    train_path = file_path + "/train/"
    train_image_files = os.listdir(train_path)
    train_images_length = len(train_image_files)
    train_images = [load_image(train_path + str(file+1) + '.tif') for file in range(train_images_length)]
    train_labels_file = open(file_path + "/train.txt", "r")
    train_labels = train_labels_file.readlines()
    train_labels_file.close()

    test_path = file_path + "/test/"
    test_image_files = os.listdir(test_path)
    test_images_length = len(test_image_files)
    test_images = [load_image(test_path + str(file+1) + '.tif') for file in range(test_images_length)]
    test_labels_file = open(file_path + "/test.txt", "r")
    test_labels = test_labels_file.readlines()
    test_labels_file.close()

    for i in range(len(train_images)):
        train_images[i] = train_images[i].astype('float32') / 255.0

    for i in range(len(test_images)):
       test_images[i] = test_images[i].astype('float32')  / 255.0
       
    train_images = np.expand_dims(train_images, axis=-1)
    train_labels = np.array(train_labels)
    train_labels = train_labels.astype(np.float32)
    train_dataset_raw = tf.data.Dataset.from_tensor_slices((train_images, train_labels))
    
    test_images = np.expand_dims(test_images, axis=-1)
    test_labels = np.array(test_labels)
    test_labels = test_labels.astype(np.float32)
    test_dataset_raw = tf.data.Dataset.from_tensor_slices((test_images, test_labels))

    for batch_size in batch_sizes:
       train_dataset = train_dataset_raw.shuffle(train_images_length).batch(batch_size)
       test_dataset = test_dataset_raw.batch(batch_size)
        
       for lr in lrs:
         for dense_layer in dense_layers:
           for layer_size in layer_sizes:
             for conv_layer in conv_layers:
                    
                name = "-{}c-{}n-{}d-{}b-{}l".format(conv_layer,layer_size,dense_layer,batch_size,lr)

                loss_object = tf.keras.losses.SparseCategoricalCrossentropy()
                optimizer = tf.keras.optimizers.Adam(lr)

                # Define our metrics
                train_loss = tf.keras.metrics.Mean('train_loss', dtype=tf.float32)
                train_accuracy = tf.keras.metrics.SparseCategoricalAccuracy('train_accuracy')
                test_loss = tf.keras.metrics.Mean('test_loss', dtype=tf.float32)
                test_accuracy = tf.keras.metrics.SparseCategoricalAccuracy('test_accuracy')

                def train_step(model, optimizer, train_images, train_labels):
                    with tf.GradientTape() as tape:
                        predictions = model(train_images, training=True)
                        loss = loss_object(train_labels, predictions)
                    grads = tape.gradient(loss, model.trainable_variables)
                    optimizer.apply_gradients(zip(grads, model.trainable_variables))

                    train_loss(loss)
                    train_accuracy(train_labels, predictions)

                def test_step(model, test_images, test_labels):
                    predictions = model(test_images)
                    loss = loss_object(test_labels, predictions)

                    test_loss(loss)
                    test_accuracy(test_labels, predictions)

                logdir = logs_base_dir + "\\" + split_type + str("{:0>5d}".format(noi*1000)) + "-" + str(no%times) + name
                train_log_dir = logdir + '\\train'
                test_log_dir = logdir + '\\test'
                train_summary_writer = tf.summary.create_file_writer(train_log_dir)
                test_summary_writer = tf.summary.create_file_writer(test_log_dir)

                # build model
                model = tf.keras.models.Sequential()
                # activation relu is better than sigmoid
                model.add(tf.keras.layers.Conv2D(layer_size, (3,3),activation=tf.nn.relu, input_shape=train_images.shape[1:]))
                model.add(tf.keras.layers.MaxPool2D(pool_size=(2,2)))

                for l in range(conv_layer-1):
                    model.add(tf.keras.layers.Conv2D(layer_size, (3,3),activation=tf.nn.relu))
                    model.add(tf.keras.layers.MaxPool2D(pool_size=(2,2)))

                model.add(tf.keras.layers.Flatten())

                for l in range(dense_layer):
                    model.add(tf.keras.layers.Dense(dense_layer_sizes[dense_layer-1-l],activation=tf.nn.relu))

                model.add(tf.keras.layers.Dense(3,activation=tf.nn.softmax))
                model.summary()

                # training for each epoch
                for epoch in range(EPOCHS):
                    for (train_images, train_labels) in train_dataset:
                        train_step(model, optimizer, train_images, train_labels)
                    with train_summary_writer.as_default():
                        tf.summary.scalar('loss', train_loss.result(), step=epoch)
                        tf.summary.scalar('accuracy', train_accuracy.result(), step=epoch)

                    for (test_images, test_labels) in test_dataset:
                        test_step(model, test_images, test_labels)
                    with test_summary_writer.as_default():
                        tf.summary.scalar('loss', test_loss.result(), step=epoch)
                        tf.summary.scalar('accuracy', test_accuracy.result(), step=epoch)

                    template = 'Epoch {}, Loss: {}, Accuracy: {}, Test Loss: {}, Test Accuracy: {}'
                    print (template.format(epoch+1,
                                           train_loss.result(), 
                                           train_accuracy.result()*100,
                                           test_loss.result(), 
                                           test_accuracy.result()*100))

                    # Reset metrics every epoch
                    train_loss.reset_states()
                    test_loss.reset_states()
                    train_accuracy.reset_states()
                    test_accuracy.reset_states()

                if noi == 20:
                    modelName = model_path + "/model_" + split_type + str("{:0>5d}".format(noi*1000)) + "-" + str(no%times) + name + ".h5"
                    model.compile(optimizer=tf.optimizers.Adam(lr),
                                  loss=tf.losses.SparseCategoricalCrossentropy(),
                                  metrics=[tf.metrics.SparseCategoricalAccuracy()])
                    model.save(modelName)

                del model

    del train_images, train_labels, test_images, test_labels, train_dataset, test_dataset
