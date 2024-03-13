import numpy as np
import os
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Flatten, Dense, Dropout
from kerastuner import BayesianOptimization

image_height = 380
image_width = 780
num_channels = 3

def load_images(folder_path):
    images = []
    labels = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".png"):
            image_path = os.path.join(folder_path, filename)
            image = tf.keras.preprocessing.image.load_img(image_path, target_size=(image_height, image_width))
            image = tf.keras.preprocessing.image.img_to_array(image)
            image = image / 255.0
            images.append(image)
            labels.append(os.path.basename(folder_path))
    return np.array(images), np.array(labels)

positives_images, positives_labels = load_images("/kaggle/input/positives-image")
negatives_images, negatives_labels = load_images("/kaggle/input/negatives-image")

positives_labels = np.ones(len(positives_images))
negatives_labels = np.zeros(len(negatives_images))

images = np.concatenate((positives_images, negatives_images), axis=0)
labels = np.concatenate((positives_labels, negatives_labels), axis=0)
indices = np.random.permutation(len(images))
images = images[indices]
labels = labels[indices]

train_images, test_images, train_labels, test_labels = train_test_split(images, labels, test_size=0.2, random_state=42)

def build_model(hp):
    model = Sequential()
    model.add(Conv2D(filters=hp.Int('conv_filters', min_value=32, max_value=128, step=32),
                     kernel_size=hp.Choice('kernel_size', values=[3, 5]),
                     activation='relu',
                     input_shape=(image_height, image_width, num_channels)))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    
    for i in range(hp.Int('num_conv_layers', 1, 3)):
        model.add(Conv2D(filters=hp.Int(f'conv_{i}_filters', min_value=32, max_value=128, step=32),
                         kernel_size=hp.Choice(f'conv_{i}_kernel_size', values=[3, 5]),
                         activation='relu'))
        model.add(MaxPooling2D(pool_size=(2, 2)))
    
    model.add(Flatten())
    
    for i in range(hp.Int('num_dense_layers', 1, 2)):
        model.add(Dense(units=hp.Int(f'dense_{i}_units', min_value=32, max_value=128, step=32),
                        activation='relu'))
        model.add(Dropout(rate=hp.Float(f'dropout_{i}_rate', min_value=0.1, max_value=0.5, step=0.1)))
    
    model.add(Dense(1, activation='sigmoid'))
    
    model.compile(optimizer=tf.keras.optimizers.Adam(hp.Choice('learning_rate', values=[1e-3, 1e-4])),
                  loss='binary_crossentropy',
                  metrics=['accuracy'])
    
    return model

tuner = BayesianOptimization(
    build_model,
    objective='val_accuracy',
    max_trials=10,
    directory='my_dir',
    project_name='cnn_tuning'
)

tuner.search(train_images, train_labels, epochs=10, validation_split=0.2)

best_model = tuner.get_best_models(num_models=1)[0]
best_model.fit(train_images, train_labels, epochs=10)
best_model.save("best_model.h5")

loaded_model = load_model("best_model.h5")
test_loss, test_accuracy = best_model.evaluate(test_images, test_labels)
print("Test Accuracy:", test_accuracy)
print("Test Loss:", test_loss)
