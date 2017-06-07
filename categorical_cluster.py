""" InterMine @ Open Genome Informatics - Similarity Project 
	-> Extension to Hierarchical Agglomerative Clustering algorithm (Single Link / Complete Link) for handling mixed numeric & categorical data type (Possiblities of multiple values in one category)          
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
			distance_matrix[row][column] = distance_mixed(first_point,second_point,numeric)


	return distance_matrix

	


""" Function to implement a Hierarchical Agglomerative Clustering Algorithm for mixed type of datasets with categorical variables being allowed to hold multiple categorical values """
def hierarchical_mixed(dataset,n_clusters,numeric):
	#Compute Distance Matrix (Complexity : O(n^2 * d)) , d=> dimension of feature
	distance_matrix = compute_distance_matrix(dataset,numeric)
	
	print distance_matrix



def create_test_data():
	#Creation of Dataset
	data = []

	#Data will have both categorical as well as numeric part
	test = [[1,2,3,4,['a','b'],[6]],[3,2,5,1,['a','b','c'],[5]]]

	#Calling the function
	hierarchical_mixed(test,2,4)
	




create_test_data()




