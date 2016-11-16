#!usr/bin/python
#Richard Wang
#Computer Simulation of Protein Folding using the 3D HP Lattice Model
#Folds protein sequence in 3D using HP Lattice Model, with Simulated Annealing

import os, sys, time
from math import *
from random import *

EMPTY = ' '

best_coords = []
final_energies = []
final_energies2 = []

protein = "PHHPPPPPPHHPPPHHHPHPPHPHHPPHPPHPPHHPPHHHHHHHPPHH"    #Protein sequence being folded



numsteps = 100000                                         #number of evaluations to perform
temperature = 600                                       #temp of simulation, 600 is good for 48 bead contact function, 1200 is good for 48 bead distance function 
starttemp = temperature
rate = float(temperature-100)/float(numsteps)		#temp decrease rate for simulated annealing

timelimit = .04						#length of simulation if running based on time, hours
threshold =  0						#energy acceptance threshold for writing if based on time, 0 for none, used to save disk space
numtrials = 3						#number of trials to run, put -1 to use time instead
starting_direction = "R"                                #R = Right, L = Left, U = Up, D = Down

energytype = 1						#Energy calculation method, -1 for distance, 1 for contact

boltzmann = 0.001987                                    #boltzmann constant(kcal/K.mol)


#Main function. Do not change anything below.

def main():
	temperature = starttemp
	element_idx = 0 # Index of the current element in protein

	coords = []             #2d array contains coords of all beads
	num_list = []                                                                                  #list that contains coord of h beads
	max_energy = 0.0                                                                                #float

	current_row = current_col = current_z = int((len(protein)*2-1)/2)                                                           #starts in middle
	num_list, coords = initialization(coords, num_list, current_row, current_col, current_z, 1)             #set up grid


	choices = [0,1,2,3,4,5]			#0,1 z axis, 2,3 x axis 4,5 y axis, one for each direction of pivot
	direction = choice(choices)

	index = randrange(1,len(protein)-1)

	new_energy =  max_energy = energy(num_list)                                  #base energy for calculations to begin
	max_energy2 = 0

	output_file(coords, 0)				                        #outputs initialization state

	current_coords = new_coords = fold(coords, num_list, direction, index)          #folds once to begin

	output_file(current_coords, 1)


	current_num_list = new_num_list =  coords_to_num_list(current_coords)
	current_energy = new_energy = energy(current_num_list)

	viable = 0								#number of possible steps
	numaccepted = 0								#number of accepted steps
	output_idx = 2								#number of times outputted to file
		
	for step in range(numsteps):
		direction = choice(choices)
		index = randrange(1,len(protein)-1)				#picks a direction
		
		new_coords = fold(current_coords, current_num_list, direction, index)		#tries a fold
		new_num_list = coords_to_num_list(new_coords)
		new_energy = energy(new_num_list)
		
		if new_num_list != current_num_list:						#if fold was possible(no overlaps)
			viable += 1
			accept = pow(e,-1*(new_energy-current_energy)/(boltzmann*temperature))	#calculate acceptance value using Metropolis algorithm
	
			rnum = random()								#gets random number between 0 and 1
			
			if rnum<accept:								#determines whether to accept or reject

				numaccepted += 1
				output_file(new_coords, output_idx)				#writes to .xyz file
				output_idx += 1
				
				global energytype
				energytype *= -1
				energy2 = energy(new_num_list)                                  #calculates comparison energy for comparing between energy functions
				energytype *= -1
				
				if (new_energy < max_energy) or ((new_energy == max_energy) and (energy2<max_energy2)):

					max_energy = new_energy
					
					energytype *= -1
					max_energy2 = energy(new_num_list)		        #calculates comparison energy with same grid, then switches back to normal
					energytype *= -1

					del best_coords[:]					#clears out previous list

					best_coords.append(new_coords)				#adds in new to list
	
				current_energy = new_energy					#adopts new coords, grid, energy, after taking it
				current_coords = new_coords
				current_num_list = new_num_list
		
		
		temperature -= rate
			

	output_file(best_coords[0], output_idx)                                                 #outputs the best fold found at the very end, for convenience when analyzing folds
	final_energies.append(str(max_energy))                                                  #ouputs best energy associated with best fold
	final_energies2.append(str(max_energy2))						#outputs conversion to alternate energy form also



#************Function Definitions**************

#initialization: places the beads into empty grid to begin folding, recursive
def initialization(coords, num_list, current_row, current_col, current_z, counter):
	if counter <= len(protein):
		 
		coords.append([(current_row), (current_col), (current_z)])

		if protein[counter-1]=='H':   							#adds coord for only 'H'
			num_list.append(current_row)
			num_list.append(current_col)
			num_list.append(current_z)

		#print (current_grid)

		if starting_direction == 'R':
			current_col += 1
		elif starting_direction == 'D':
			current_row += 1
		elif starting_direction == 'L':
			current_col -= 1
		elif starting_direction == 'U':
			current_row -= 1                
		elif starting_direction == 'F':
			current_z += 1
		elif starting_direction == 'B':
			current_z -= 1

		num_list, coords = initialization(coords, num_list, current_row, current_col, current_z, counter+1)

	return num_list, coords

	

#energy: calculates energy for current grid, no movement
def energy(num_list):
	energy = 0
	if energytype == -1:            #distance
		if len(num_list) >3:                                                 #avoids redundant calling
			for i in range ((len(num_list)/3)-1):
				for j in range (i+1,(len(num_list)/3)):            #avoids duplicates
					xdif = num_list[3*i]-num_list[3*j]                #finds coord differences
					ydif = num_list[3*i+1]-num_list[3*j+1]
	
					zdif = num_list[3*i+2]-num_list[3*j+2]
	
					energy -= 1/(sqrt(pow(xdif,2) + pow(ydif,2) + pow(zdif,2)))           #adds inverse energy using distance calculation formula
		return energy
	
	elif energytype == 1:           #contact
		if len(num_list) >3:                                                 #avoids redundant calling
			for i in range ((len(num_list)/3)-1):
				for j in range (i+1,(len(num_list)/3)):
					z = num_list[3*j+2]
					y = num_list[3*j+1]
					x = num_list[3*j]

					x2 = num_list[3*i]
					y2 = num_list[3*i+1]
					z2 = num_list[3*i+2]

					if (abs(x-x2)==1 and y == y2 and z == z2) or (x == x2 and abs(y-y2)==1 and z == z2) or (x == x2 and y == y2 and abs(z-z2)==1):
						energy -= 1
						
		return energy


#fold: performs a folding maneuever given available at bead num element_idx pointing at direction
def fold(coords, num_list, direction, element_idx):
	#direction is either 0,1,2,3,4,5. Corresponding directions listed above.
	#element_idx is pivot point
	
	translated_coords = []
	for i in range(len(coords)):                                            #copies the nested array
		translated_coords.append([])
		translated_coords[i].append(coords[i][0])
		translated_coords[i].append(coords[i][1])
		translated_coords[i].append(coords[i][2])

	del translated_coords [0: element_idx]                                  #removes beginning of coords, keeps pivot point coords
	trans_x = translated_coords[0][0]
	trans_y = translated_coords[0][1]
	trans_z = translated_coords[0][2]
	
	for i in range (len(translated_coords)):                                #shifts all coords so that origin is at 0,0 and others are relative
		translated_coords[i][0] -= trans_x
		translated_coords[i][1] -= trans_y
		translated_coords[i][2] -= trans_z	

	rotated_coords = []

	for i in range (len(translated_coords)):                                            #copies the nested array
		rotated_coords.append([])
		rotated_coords[i].append(translated_coords[i][0])
		rotated_coords[i].append(translated_coords[i][1])
		rotated_coords[i].append(translated_coords[i][2])

		#************Rotation Calculations****************

		if direction == 0:		#z axis
			angle = pi/2
			rotated_coords[i][0] = int(round(translated_coords[i][0] * cos(angle) - translated_coords[i][1] * sin(angle),0))        #rotation formula
			rotated_coords[i][1] = int(round(translated_coords[i][0] * sin(angle) + translated_coords[i][1] * cos(angle),0))
		
		elif direction == 1:
			angle = -pi/2
			rotated_coords[i][0] = int(round(translated_coords[i][0] * cos(angle) - translated_coords[i][1] * sin(angle),0))        #rotation formula
			rotated_coords[i][1] = int(round(translated_coords[i][0] * sin(angle) + translated_coords[i][1] * cos(angle),0))

		elif direction == 2:		#x axis
			angle = pi/2
			rotated_coords[i][1] = int(round(translated_coords[i][1] * cos(angle) - translated_coords[i][2] * sin(angle),0))        #rotation formula
			rotated_coords[i][2] = int(round(translated_coords[i][1] * sin(angle) + translated_coords[i][2] * cos(angle),0))

		elif direction == 3:
			angle = -pi/2
			rotated_coords[i][1] = int(round(translated_coords[i][1] * cos(angle) - translated_coords[i][2] * sin(angle),0))        #rotation formula
			rotated_coords[i][2] = int(round(translated_coords[i][1] * sin(angle) + translated_coords[i][2] * cos(angle),0))

		elif direction == 4:		#y axis
			angle = pi/2
			rotated_coords[i][2] = int(round(translated_coords[i][2] * cos(angle) - translated_coords[i][0] * sin(angle),0))        #rotation formula
			rotated_coords[i][0] = int(round(translated_coords[i][2] * sin(angle) + translated_coords[i][0] * cos(angle),0))

		elif direction == 5:
			angle = -pi/2
			rotated_coords[i][2] = int(round(translated_coords[i][2] * cos(angle) - translated_coords[i][0] * sin(angle),0))        #rotation formula
			rotated_coords[i][0] = int(round(translated_coords[i][2] * sin(angle) + translated_coords[i][0] * cos(angle),0))

		
		#************End Rotation Calculations***************

		rotated_coords[i][0] += trans_x                                 #shifts coords back to original relative position
		rotated_coords[i][1] += trans_y
		rotated_coords[i][2] += trans_z

	final_coords = []

	for i in range (len(coords)):                                            #copies the nested array once again to append on rotated coords
		if i < element_idx:
			final_coords.append([])
			final_coords[i].append(coords[i][0])
			final_coords[i].append(coords[i][1])
			final_coords[i].append(coords[i][2])
		else:
			final_coords.append([])					#copies the new part over to the old beginning
			final_coords[i].append(rotated_coords[i-element_idx][0])
			final_coords[i].append(rotated_coords[i-element_idx][1])
			final_coords[i].append(rotated_coords[i-element_idx][2])
			
	overlap = 0

	for i in range (0,len(coords)):                                         #checks for overlap before accepting move, if overlap it will return the original without fold
		for j in range (element_idx+1,len(final_coords)):
			#print (coords[i][0], final_coords[j][0], coords[i][1], final_coords[j][1], i, j)
			if coords[i][0] == final_coords[j][0] and coords[i][1] == final_coords[j][1] and coords[i][2] == final_coords[j][2]:    #detects overlap to fulfill self-avoiding walk
				return coords                   #If overlap detected, returns original pre-fold conformation
	return final_coords                                     #If no overlap detected, returns the post-fold conformation
		
#coords_to_num_list: converts the coords data into num_list data
def coords_to_num_list(coords):
	num_list = []
	for i in range(len(protein)):
		if protein[i] == "H":
			num_list.append(coords[i][0])                           #adds H beads only to num_list, no P beads
			num_list.append(coords[i][1])
			num_list.append(coords[i][2])
	return num_list

#output_file: writes into file 'output.xyz' visualizable via vmd
def output_file(coords,output_idx):
	
	trans_x = coords[0][0]					#retrieves first bead coordinates to center protein around (0,0,0)
	trans_y = coords[0][1]
	trans_z = coords[0][2]

	f.write(str(len(protein))+"\n")
	f.write("Frame "+str(output_idx)+"\n")
	for i in range(len(protein)):
		f.write(protein[i]+"\t"+str(float(coords[i][1]-trans_y))+"\t"+str(-1 * float(coords[i][0]-trans_x))+"\t"+str(-1 * float(coords[i][2]-trans_z))+"\n")			#writes to .xyz file in specific format



#******* Calls main function *******

#Main program
start_time = time.time()		#timer keeps track of execution time

if numtrials != -1:                     #If trial based data collection
	for trial in range(numtrials):
		f = open ("output_%s.xyz"%(trial),"w+") #opens file to write in
		f.flush()				#erases previous contents if any
		
		main()					#main program run, writes coordinate data into file f
		final_energies[trial] = str(trial)+"  "+final_energies[trial]+"  "+final_energies2[trial]	

		f.close()                               #closes file
	

	f = open ("data.txt","w+")					#general data output file
	f.flush()                                                       #erases previous contents if any

	s = ""
	for i in range (len(final_energies)):				#outputs final energies of each trial to file if running in batches
		s += final_energies[i] + "\n"

	f.write(s)
	f.write(str(time.time()-start_time))				#outputs final runtime to file
	f.close()

else:					#timed data collection mode, will filter out results to only below a certain threshold specified to save storage space
	trial = 0                       #indexes for output, number of trials passed threshold
	totaltrials = 0                 #indexes for output, number of total trials

	d = open ("data.txt","w+")
	d.close()

	while ((time.time()-start_time)/60.0/60.0)<timelimit:
		f = open("output_%s.xyz"%(trial),"w+")
		f.truncate()
		main()
		f.close()
	
		#print final_energies[totaltrials]
		#print int(round(float(final_energies[totaltrials]),0)), threshold
		if int(round(float(final_energies[totaltrials]),0)) < threshold:
			final_energies[totaltrials] = str(totaltrials) + "  " + str(trial)+"  "+final_energies[totaltrials]+"  "+final_energies2[totaltrials]
			s = ""
			s += final_energies[totaltrials] + "\n"
			d = open ("data.txt","a")
			d.write(s)                                      #writes data into file if past threshold

			d.close()		
		else:
			os.remove("output_%s.xyz"%(trial))              #otherwise removes the file to save space
			trial -= 1		
		trial += 1
		totaltrials += 1
	d = open ("data.txt","a")
	d.write(str(totaltrials))                                       #outputs total number of trials at the end
	d.close()

#end
