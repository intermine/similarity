""" InterMine @ Open Genome Informatics - Similarity Project 
	-> Extension to Hierarchical Agglomerative Clustering algorithm  for handling mixed numeric & categorical data type (Possiblities of multiple values in one category)          
	-> The distance metric will be divided into two parts : For Numeric => Euclidean, For categories => Jaccard Coefficient   
	-> Based on http://edu.cs.uni-magdeburg.de/EC/lehre/sommersemester-2013/wissenschaftliches-schreiben-in-der-informatik/publikationen-fuer-studentische-vortraege/kMeansMixedCatNum.pdf """



#Libraries
from __future__ import division
import json
import pandas as pd 
import numpy as np
import string
import re
import sys
import matplotlib as plt 
import random
import math


#Function to visualize the clustering as a dendogram - To be completed via formation of linkage matrix
def dendogram():
	return 1

#Function to compute intersection of two lists
def intersect(list1,list2):
	return list(set(list1) & set(list2))

#Function to compute union of two lists
def union(list1,list2):
	return list(set(list1) | set(list2))

""" _Function_ name : distance_mixed 
	 @Parameters : 
		 1. datapoint1 => First feature vector
		 2. datapoint2 => Second feature vector
		 3. numeric => Number of numeric features in the vector 

	 @Return : Distance b/w two mixed data points                                """

def distance_mixed(datapoint1,datapoint2,numeric):
	#Numeric part of Two Points
	point1_numeric = datapoint1[:numeric]
	point2_numeric = datapoint2[:numeric]

	#Conversion into Numpy Array	
	point1_numeric = np.array(point1_numeric)
	point2_numeric = np.array(point2_numeric)

	#Distance b/w numeric features
	euclidean_distance = np.linalg.norm(point1_numeric - point2_numeric)

	#Categorical Part of Two Points
	point1_categorical = datapoint1[numeric:]
	point2_categorical = datapoint2[numeric:]

	jaccard_distance = 0

	#Computing Jaccard Distance
	for feature1,feature2 in zip(point1_categorical,point2_categorical):
		#Intersection
		intersection_feature = intersect(feature1,feature2)
		#Union
		union_feature = union(feature1,feature2)

		jaccard_distance += 1 - len(intersection_feature)/len(union_feature)


	return (euclidean_distance + jaccard_distance)


#Function to compute distances b/w every pair of points
def compute_distance_matrix(dataset,numeric):
	#Size of the distance matrix
	n = len(dataset)

	#Initializing the matrix with zero
	distance_matrix = np.zeros(shape = (n,n))
	
	#Computation of the distance matrix
	for first_point in dataset:
		row = dataset.index(first_point)
		for second_point in dataset:
			column = dataset.index(second_point)
			if row == column:
				#As same elements are in the same cluster itself -- To avoid minimum distance computation
				distance_matrix[row][column] = float("inf")
			else:
				distance_matrix[row][column] = distance_mixed(first_point,second_point,numeric)


	return distance_matrix

	
#Function to find the minimum distance b/w two datapoints
def find_distance_min(distance_matrix):
	#Minimum value in the matrix
	minimum = np.min(distance_matrix)
	#Indices of the Minimum Value - List of indices holding minimum value
	position = np.argwhere(distance_matrix == minimum)

	#Extract the first element of the list
	position = position[0]
	#Row
	row = position[0]
	#Column
	column = position[1]
	#Reset the position in the matrix
	distance_matrix[row][column] = float("inf")
	distance_matrix[column][row] = float("inf")

	return distance_matrix,row ,column



""" _Function_ name : hierarchical_mixed 
	   Objective : Unsupervised Hierarchical Clustering Algorithm which is able to cluster data based on mixed types 
	   @Parameters : 
		   dataset : List of feature vectors
		   n_clusters : Number of desired clusters 
		   numeric : Number of numeric features in the feature vector

	   @Return : Dictionary with feature vector number as keys and their corresponding cluster as values                    """

def hierarchical_mixed(dataset,n_clusters,numeric):
	#Compute Initial Distance Matrix (Complexity : O(n^2 * d)) , d=> dimension of feature
	distance_matrix = compute_distance_matrix(dataset,numeric)      

	#Initial Condition : All the points are individual clusters
	number_clusters = len(dataset)

	#Dictionary for cluster assignment
	clusters = {}

	cluster_no = 0 

	#dataset = list(set(dataset))

	print number_clusters
	#Initialization of Clusters
	
	for feature_vector in dataset:
		#clusters[dataset.index(feature_vector)] = cluster_no
		clusters[cluster_no] = cluster_no
		cluster_no += 1

	print clusters

	#Merging Process
	while number_clusters > n_clusters:
		#Find the Minimum Distance and Position of the two points to be merged
		#print 1
		distance_matrix, first_point,second_point = find_distance_min(distance_matrix)
		#Temporary points for distance matrix manipulation
		temp_point1 = clusters[first_point]
		temp_point2 = clusters[second_point]

		#Manipulation of distance matrix for resetting distances b/w points in same cluster as infinity
		for datapoint1 in clusters:
			if clusters[datapoint1] == temp_point1:
				for datapoint2 in clusters:
					if clusters[datapoint2] == temp_point2:
						distance_matrix[datapoint1][datapoint2] = float("inf")
						distance_matrix[datapoint2][datapoint1] = float("inf")


        #Merging two clusters together
		for datapoint in clusters:
			if clusters[datapoint] == temp_point2:
				#Reasignment
				clusters[datapoint] = clusters[first_point]

		#Manipulation of the Distance Matrix


		#clusters[second_point] = clusters[first_point]
		number_clusters -=1	

    #Returning a dictionary corresponding to node ID and their cluster number
	return clusters
	
	


def create_test_data():
	#Creation of Dataset
	data = []

	#Data will have both categorical as well as numeric part
	test = [[1,2,3,4,['a','b'],[6]],[3,2,5,1,['a','b','c'],[5]],[1,2,3,4,['a','b','c'],[6]],[3,1,5,2,['a','b'],[5]]]

	#Calling the function
	hierarchical_mixed(test,2,4)
	




#create_test_data()




