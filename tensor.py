""" Denoising Autoencoder for Dimensionality Reduction """

#Libraries
import tensorflow as tf 
import numpy as np 
import math


#Function to calculate loss
def calculate_loss(predicted,actual):
	#Cross Entropy error
	cross_entropy = -tf.reduce_sum(actual*tf.log(predicted))

	return cross_entropy


#AutoEncoder Definition
def AutoEncoder(input):
	#Input Sample Random 
	#input = np.random.rand(1000,20)

	#Conversion into float32 for homogeneous matrix multiplication
	input = input.astype(np.float32,copy=False)

	#Adding noise to the Sample -- Perturbed Input
	noisy_input = input + 0.2 * np.random.random_sample((input.shape))

	#Scale your input data to [0,1]
	scaled_input_data = np.divide((noisy_input - noisy_input.min()),(noisy_input.max() - noisy_input.min()))

	#For mapping output as input while training
	output = input 

	#Scale your output data to [0,1]
	scaled_output_data = np.divide((output - output.min()),(output.max() - output.min()))

	input = scaled_input_data

	output = scaled_output_data

	#Initialize number of neurons in the hidden layer
	n_hidden = 2 

	#Number of Samples
	n_samples = input.shape[0]

	#Number of Input Features
	n_input = input.shape[1]

	#Input which will be fed
	x = tf.placeholder(tf.float32,[None,n_input])

	#Initializing the weights for the hidden layer
	W_hidden = tf.Variable(tf.random_uniform((n_input,n_hidden)))

	#Initializing the biases for the hidden layer
	b_hidden = tf.Variable(tf.zeros([n_hidden]))

	#Hidden Layer values
	h = tf.nn.tanh(tf.matmul(x,W_hidden) + b_hidden)

	#Second layer Initialization
	W_output = tf.transpose(W_hidden)

	#Output Layer Biases
	b_output = tf.Variable(tf.zeros([n_input]))

	#Output values 
	predicted = tf.nn.tanh(tf.matmul(h,W_output) + b_output)

	#Placeholder for output
	actual = tf.placeholder(tf.float32,[None,n_input])

	#Find the loss value
	loss = calculate_loss(predicted,actual)

	#Mean sq loss
	#loss = tf.reduce_mean(tf.square(actual - predicted))

	#Training Step
	training = tf.train.GradientDescentOptimizer(0.05).minimize(loss)

	#Initialize all variables
	init = tf.initialize_all_variables()

	with tf.Session() as sess:
		#Initialize the variables in this Session
		sess.run(init)

		#Number of Rounds
		n_rounds = 1000

		#Selection of batch size of 50
		batch_size = min(n_samples,100)

		for i in range(0,n_rounds):
			#Sample Selection -- From the total sample - samples corresponding to batch size is selected
			sample = np.random.randint(n_samples,size=batch_size)
			
			#Input data corresponding to batch size
			input_x = input[sample][:]

			#Output data corresponding to batch size
			output_y = output[sample][:]

			#Run the training
			sess.run(training, feed_dict={x:input_x,actual:output_y})

			if i%50 == 0:
				print sess.run(loss,feed_dict={x:input_x,actual:output_y})

		weights =  sess.run(W_hidden)
		biases = sess.run(b_hidden)


	return weights, biases	


		





#AutoEncoder(np.random.rand(1000,20))






