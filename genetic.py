"""  Implementation of Genetic Algorithm for Optimization for Clustering    
     To be integrated with K-means for finding initial seeds to get better results               """

from __future__ import division
import math
from random import randint, random


#Function to create an individual in a population
def create_individual(count,minimum,maximum):
	#Create random individual in a Population
	individual = [randint(minimum,maximum) for i in range(count)]

	return individual


#Function to define the fitness of the individual -- Fitness Function
def fitness(individual,target):
	#Sum of elements in the list
	individual_sum = sum(individual)
	#Distance b/w target and individual -- Fitness
	difference = abs(individual_sum - target)

	return difference

#Function to create the total population
def create_population(number_individuals,count,minimum,maximum):
	#Generate Population
	population = [create_individual(count,minimum,maximum) for individual in range(number_individuals)]

	return population


#Function to evolve the population
def evolve(population,number_individuals,count,minimum,maximum,top,mutation,random_select,target):

	#Find the fitness of each individual
	fitness_population = [(fitness(individual,target),individual) for individual in population]

	#Sort the population
	sorted_population = sorted(fitness_population)

	#Extracting the lists from the tuples to store as a parent
	parents = [x[1] for x in sorted_population]

	#Select the best fit individuals ~ 20% of the population is selected
	best_fit = parents[:int(number_individuals*top)]

	#Randomly select individuals from the unfit portion to add to the parents
	#This ensures Genetic Diversity and also ensures that the solution does not get stuck in a local maximum
	for individual in parents[int(number_individuals*top):]:
		#Randomly select an individual
		if random_select > random():
			#Choose the individual
			parents.append(individual)


	#Produce some mutation into the parents population now
	#Why Mutation? : Mutation Divergence Phenomenon, hence should be kept low
	#Mutation induces genetic diversity and hence tries to pull out of local minimas/maximas and extract better ones
	for individual in parents:
		#Keep mutation probability low -- as mutation is divergence phenomenon, not convergence
		if mutation > random():
			#Determine the position to mutate
			position_to_mutate = randint(0,len(individual)-1)

			#Replace the values in that position with some other value
			individual[position_to_mutate] = randint(min(individual),max(individual))


	#Creation of progenies / children 
	#This step is essential for convergence and tries to find a local maxima / minima 
	children_number = number_individuals - len(parents)

	children = []
    
    #Create children by merging from the parents
	while len(children) < children_number:
		#Selecting one male randomly
		male_index = randint(0,len(parents)-1)

		#Selecting one female randomly
		female_index = randint(0,len(parents)-1)

        #Ensure male and female both aren't the same
		if male_index != female_index:
			#Male
			male = parents[male_index]
			#Female
			female = parents[female_index]

			division = len(male) / 2
            
            #New Children
			new_children = male[:division] + female[division:]

			#Append to children list
			children.append(new_children)

    
    #Appending to main parents as iterables
	parents.extend(children)

	return parents


#Function to run cycles of the Evolutionary Algorithm
def main(number_individuals,count,minimum,maximum):
	#Initial Population Generation
	population = create_population(number_individuals,count,minimum,maximum)

	#Run iterations for the EA
	for i in range(0,9):
		get_parents = evolve(population,100,5,1,100,0.2,0.01,0.05,200)
		population = get_parents


	fitness_population = [(fitness(individual,200),individual) for individual in population]

	sorted_fitness_population =  sorted(fitness_population)

	solution = sorted_fitness_population[0]

	return solution




answer = main(100,5,1,100)
































	
