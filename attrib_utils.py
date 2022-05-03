from __future__ import absolute_import
import numpy as np
from keras import backend as K
import matplotlib.pyplot as plt
import tensorflow as tf
import scipy.ndimage as ndimage

def model_modifier_function(cloned_model):
    cloned_model.layers[1].activation = tf.keras.activations.linear
def macro_f1(y, y_hat, thresh=0.5):
    """Compute the macro F1-score on a batch of observations (average F1 across labels)

    Args:
        y (int32 Tensor): labels array of shape (BATCH_SIZE, N_LABELS)
        y_hat (float32 Tensor): probability matrix from forward propagation of shape (BATCH_SIZE, N_LABELS)
        thresh: probability value above which we predict positive

    Returns:
        macro_f1 (scalar Tensor): value of macro F1 for the batch
    """
    y_pred = tf.cast(tf.greater(y_hat, thresh), tf.float32)
    tp = tf.cast(tf.math.count_nonzero(y_pred * y, axis=0), tf.float32)
    fp = tf.cast(tf.math.count_nonzero(y_pred * (1 - y), axis=0), tf.float32)
    fn = tf.cast(tf.math.count_nonzero((1 - y_pred) * y, axis=0), tf.float32)
    f1 = 2 * tp / (2 * tp + fn + fp + 1e-16)
    macro_f1 = tf.reduce_mean(f1)
    return macro_f1
def macro_soft_f1(y, y_hat):
    """Compute the macro soft F1-score as a cost.
    Average (1 - soft-F1) across all labels.
    Use probability values instead of binary predictions.

    Args:
        y (int32 Tensor): targets array of shape (BATCH_SIZE, N_LABELS)
        y_hat (float32 Tensor): probability matrix of shape (BATCH_SIZE, N_LABELS)

    Returns:
        cost (scalar Tensor): value of the cost function for the batch
    """

    y = tf.cast(y, tf.float32)
    y_hat = tf.cast(y_hat, tf.float32)
    tp = tf.reduce_sum(y_hat * y, axis=0)
    fp = tf.reduce_sum(y_hat * (1 - y), axis=0)
    fn = tf.reduce_sum((1 - y_hat) * y, axis=0)
    soft_f1 = 2 * tp / (2 * tp + fn + fp + 1e-16)
    cost = 1 - soft_f1  # reduce 1 - soft-f1 in order to increase soft-f1
    macro_cost = tf.reduce_mean(cost)  # average on all labels

    return macro_cost
def normalize(array, min_value=0., max_value=1.):
    """Normalizes the numpy array to (min_value, max_value)
    Args:
        array: The numpy array
        min_value: The min value in normalized array (Default value = 0)
        max_value: The max value in normalized array (Default value = 1)
    Returns:
        The array normalized to range between (min_value, max_value)
    """
    arr_min = np.min(array)
    arr_max = np.max(array)
    normalized = (array - arr_min) / (arr_max - arr_min + K.epsilon())
    return (max_value - min_value) * normalized + min_value
def get_img_shape(img):
    """Returns image shape in a backend agnostic manner.
    Args:
        img: An image tensor of shape: `(channels, image_dims...)` if data_format='channels_first' or
            `(image_dims..., channels)` if data_format='channels_last'.
    Returns:
        Tuple containing image shape information in `(samples, channels, image_dims...)` order.
    """
    if isinstance(img, np.ndarray):
        shape = img.shape
    else:
        shape = K.int_shape(img)

    if K.image_data_format() == 'channels_last':
        shape = list(shape)
        shape.insert(1, shape[-1])
        shape = tuple(shape[:-1])
    return shape
def deprocess_input(input_array, input_range=(0, 255)):
    """Utility function to scale the `input_array` to `input_range` throwing away high frequency artifacts.
    Args:
        input_array: An N-dim numpy array.
        input_range: Specifies the input range as a `(min, max)` tuple to rescale the `input_array`.
    Returns:
        The rescaled `input_array`.
    """
    # normalize tensor: center on 0., ensure std is 0.1
    input_array = input_array.copy()
    input_array -= input_array.mean()
    input_array /= (input_array.std() + K.epsilon())
    input_array *= 0.1

    # clip to [0, 1]
    input_array += 0.5
    input_array = np.clip(input_array, 0, 1)

    # Convert to `input_range`
    return (input_range[1] - input_range[0]) * input_array + input_range[0]
def random_array(shape, mean=128., std=20.):
    """Creates a uniformly distributed random array with the given `mean` and `std`.
    Args:
        shape: The desired shape
        mean: The desired mean (Default value = 128)
        std: The desired std (Default value = 20)
    Returns: Random numpy array of given `shape` uniformly distributed with desired `mean` and `std`.
    """
    x = np.random.random(shape)
    # normalize around mean=0, std=1
    x = (x - np.mean(x)) / (np.std(x) + K.epsilon())
    # and then around the desired mean/std
    x = (x * std) + mean
    return x

    """Ensures that the value is a list. If it is not a list, it creates a new list with `value` as an item.
    """
    if not isinstance(value, list):
        value = [value]
    return value


# saliency maps
def saliency_map_dense(model, test_field):

    images = tf.Variable(test_field, dtype=float)
    #print(images.shape)
    #tf.reshape(images,[None, 30, 30, 30, 3])
    with tf.GradientTape() as tape:
        pred = model(images, training=False)
        #np.array(pred)
        class_idxs_sorted = np.argsort(pred.numpy().flatten())[::-1]
        loss = pred[0][class_idxs_sorted[0]]    
        grads = tape.gradient(loss, images)

    dgrad_abs = tf.math.abs(grads)
    #print(dgrad_abs.shape)

    #dgrad_max_ = np.max(dgrad_abs, axis=3)#[0]

    #print(dgrad_max_.shape)

    ## normalize to range between 0 and 1
    arr_min, arr_max  = np.min(dgrad_abs), np.max(dgrad_abs)
    grad_eval = (dgrad_abs - arr_min) / (arr_max - arr_min + 1e-18)
    #print(grad_eval)
    #print(grad_eval.shape)
    #plt.hist(grad_eval.flatten())
    #plt.show()
    return np.around(grad_eval/10, decimals=2)*10


# integrated gradients
def get_gradients(img_input, model, top_pred_idx):
    """Computes the gradients of outputs w.r.t input image.

    Args:
        img_input: 4D image tensor
        top_pred_idx: Predicted label for the input image

    Returns:
        Gradients of the predictions w.r.t img_input
    """
    images = tf.cast(img_input, tf.float32)

    with tf.GradientTape() as tape:
        tape.watch(images)
        preds = model(images)
        top_class = preds[:, top_pred_idx]

    grads = tape.gradient(top_class, images)
    return grads
def get_integrated_gradients(img_input, top_pred_idx, model, baseline=None, num_steps=50):
    """Computes Integrated Gradients for a predicted label.

    Args:
        img_input (ndarray): Original image
        top_pred_idx: Predicted label for the input image
        baseline (ndarray): The baseline image to start with for iterpolation
        num_steps: Number of interpolation steps between the baseline
            and the input used in the computation of integrated gradients. These
            steps along determine the integral approximation error. By default,
            num_steps is set to 50.

    Returns:
        Integrated gradients w.r.t input image
    """
    img_size = np.shape(img_input)
    # If baseline is not provided, start with a black image
    # having same size as the input image.
    if baseline is None:
        baseline = np.zeros(img_size).astype(np.float32)
    else:
        baseline = baseline.astype(np.float32)

    # 1. Do interpolation.
    img_input = img_input.astype(np.float32)
    interpolated_image = [
        baseline + (step / num_steps) * (img_input - baseline)
        for step in range(num_steps + 1)
    ]
    interpolated_image = np.array(interpolated_image).astype(np.float32)

    # 2. Preprocess the interpolated images
    #interpolated_image = xception.preprocess_input(interpolated_image)

    # 3. Get the gradients
    grads = []
    for i, img in enumerate(interpolated_image):
        img = tf.expand_dims(img, axis=0)
        grad = get_gradients(img, model, top_pred_idx=top_pred_idx)
        grads.append(grad[0])
    grads = tf.convert_to_tensor(grads, dtype=tf.float32)

    # 4. Approximate the integral using the trapezoidal rule
    grads = (grads[:-1] + grads[1:]) / 2.0
    avg_grads = tf.reduce_mean(grads, axis=0)

    # 5. Calculate integrated gradients and return
    integrated_grads = (img_input - baseline) * avg_grads
    return integrated_grads
def random_baseline_integrated_gradients(
    img_input, model, top_pred_idx, num_steps=50, num_runs=2
):
    """Generates a number of random baseline images.

    Args:
        img_input (ndarray): 3D image
        top_pred_idx: Predicted label for the input image
        num_steps: Number of interpolation steps between the baseline
            and the input used in the computation of integrated gradients. These
            steps along determine the integral approximation error. By default,
            num_steps is set to 50.
        num_runs: number of baseline images to generate

    Returns:
        Averaged integrated gradients for `num_runs` baseline images
    """
    # 1. List to keep track of Integrated Gradients (IG) for all the images
    integrated_grads = []
    img_size = np.shape(img_input)
    # 2. Get the integrated gradients for all the baselines
    for run in range(num_runs):
        #baseline = np.random.random(img_size) * 255
        baseline = random_array(img_size, mean = 0, std = 1)
        igrads = get_integrated_gradients(
            img_input=img_input,
            model = model, 
            top_pred_idx=top_pred_idx,
            baseline=baseline,
            num_steps=num_steps,
        )
        integrated_grads.append(igrads)

    # 3. Return the average integrated gradients for the image
    integrated_grads = tf.convert_to_tensor(integrated_grads)
    return tf.reduce_mean(integrated_grads, axis=0)



# grad-cam
def make_gradcam_heatmap(img_array, model, last_conv_layer_name, pred_index=None):
    # First, we create a model that maps the input image to the activations
    # of the last conv layer as well as the output predictions
    grad_model = tf.keras.models.Model(
        [model.inputs], [model.get_layer(last_conv_layer_name).output, model.output]
    )

    # Then, we compute the gradient of the top predicted class for our input image
    # with respect to the activations of the last conv layer
    with tf.GradientTape() as tape:
        last_conv_layer_output, preds = grad_model(img_array)
        if pred_index is None:
            pred_index = tf.argmax(preds[0])
        class_channel = preds[:, pred_index]

    # This is the gradient of the output neuron (top predicted or chosen)
    # with regard to the output feature map of the last conv layer
    grads = tape.gradient(class_channel, last_conv_layer_output)

    # This is a vector where each entry is the mean intensity of the gradient
    # over a specific feature map channel
    pooled_grads = tf.reduce_mean(grads, axis=(0, 1, 2))

    # We multiply each channel in the feature map array
    # by "how important this channel is" with regard to the top predicted class
    # then sum all the channels to obtain the heatmap class activation
    last_conv_layer_output = last_conv_layer_output[0]
    heatmap = last_conv_layer_output @ pooled_grads[..., tf.newaxis]
    heatmap = tf.squeeze(heatmap)

    # For visualization purpose, we will also normalize the heatmap between 0 & 1
    heatmap = tf.maximum(heatmap, 0) / tf.math.reduce_max(heatmap)
    return heatmap.numpy()