#This script implements the velocity-verlet algorithm
#The argon is used here and the mass =39.948
# sigma=0.3345  epsilon=0.0661

import numpy as np
import sys
import os

#checked: correct
def lj(r):
	sigma=0.3345  
	epsilon=0.0661
	return 4*epsilon*((sigma/r)**12-(sigma/r)**6)

#checked:correct
def V_lj(coord):
	distances=get_distances(coord)
	sigma=0.3345  
	epsilon=0.0661
	cum=0
	row,col=np.shape(distances)
	for i in range(row):
		for j in range(i+1,col):
			r=distances[i,j]
			if r>0:
				cum=cum+4*epsilon*((sigma/r)**12-(sigma/r)**6)
			else:
				print("bad contact detected.")
	return cum

#checked:correct
def read_Natoms(initial_file):
	content=open(initial_file,"r")
	for i,lines in enumerate(content):
		if i==0:
			#print(lines)
			Natoms=int(lines)
			mass=np.zeros(Natoms)
			coord=np.zeros((Natoms,3))
			k=0
		elif i==1:
			pass
		else:
			#the element_name to mass need to implemented in the future.
			mass_name,coord[k,0],coord[k,1],coord[k,2]=lines.split()
			mass[k]=39.948
			k+=1
	return Natoms,mass,coord

#checked:correct
def get_distances(coord):
	row,col=np.shape(coord)
	matrix=np.zeros((row,row))
	for i in range(row):
		coord_i=np.zeros((row,3))
		coord_i[:,0]=coord[i,0]
		coord_i[:,1]=coord[i,1]
		coord_i[:,2]=coord[i,2]
		matrix[i:,]=np.sqrt( np.sum( (coord-coord_i)*(coord-coord_i) ,axis=1) )
	return matrix

def E_tot(coord,mass,velocity):
	vel=np.array(velocity)
	print("dongneng:",np.sum(  np.sum(vel*vel,axis=1)*mass*0.5  )) 
	return np.sum(np.sum(vel*vel,axis=1)*mass*0.5)+V_lj(coord)

#checked:correct
def get_acceleration(coord, mass):
	delta=0.001
	row,col=np.shape(coord)
	acceleration=np.zeros_like(coord)
	delta_x=np.zeros_like(coord)
	delta_y=np.zeros_like(coord)
	delta_z=np.zeros_like(coord)
	for i in range(row):
		delta_x[i,0]=delta
		delta_y[i,1]=delta
		delta_z[i,2]=delta
		acceleration[i,0]=(V_lj(coord+delta_x)-V_lj(coord-delta_x))/(2*delta)/mass[i]/(-1)
		acceleration[i,1]=(V_lj(coord+delta_y)-V_lj(coord-delta_y))/(2*delta)/mass[i]/(-1)
		acceleration[i,2]=(V_lj(coord+delta_z)-V_lj(coord-delta_z))/(2*delta)/mass[i]/(-1)
		delta_x[i,0]=0
		delta_y[i,1]=0
		delta_z[i,2]=0
	return acceleration

def verlet(coord,  mass,  velocity,Nsteps=1000,dt=0.1):
	output=open("movie.xyz","w")
	for i in range(Nsteps):
		print(i)
		#print("step "+str(i)+"   "+str(E_tot(coord,mass,velocity)))
		#print(coord)
		output.write(str(len(mass)) + "\n Frame:"+str(i)+"\n")
		for j in range(len(mass)):
			output.write("Ar\t"+str(coord[j,0]) +"\t"+str(coord[j,1]) +"\t"+str(coord[j,2])+"\n")
		vel=np.array(velocity) #convert the list to numpy.ndarray
		new_coord=coord+vel*dt+get_acceleration(coord, mass)*dt*dt/2.0
		new_velocity=vel +(get_acceleration(coord, mass)+get_acceleration(new_coord, mass))*dt/2.0
		coord=new_coord
		velocity=new_velocity
	output.close()
	return new_coord,new_velocity

def leap_frog(coord,  mass,  velocity,Nsteps=1000,dt=0.1):
	output=open("movie1.xyz","w")
	for i in range(Nsteps):
		print(i)
		output.write(str(len(mass)) + "\n Frame:"+str(i)+"\n")
		for j in range(len(mass)):
			output.write("Ar\t"+str(coord[j,0]) +"\t"+str(coord[j,1]) +"\t"+str(coord[j,2])+"\n")
		vel=np.array(velocity) #convert the list to numpy.ndarray #the initial velocity,time t-delatt
		new_velocity=vel+get_acceleration(coord, mass)*dt   #the initial velocity,time t+delatt
		new_coord=coord+new_velocity*dt
		#velocity at time t
		velocity_t=0.5*(vel+new_velocity)
		#calculate the E_tot at time t
		#print("step "+str(i)+"   "+str(E_tot(coord,mass,velocity)))
		#print(coord)

		coord=new_coord
		velocity=new_velocity
	output.close()
	return new_coord,new_velocity

def main():
	initial_file=sys.argv[1]
	Natoms,mass,coord=read_Natoms(initial_file)
	distances=get_distances(coord)
    # print("reading initial file", initial_file)
	#random velocity distrition (0,10)
	#velocity = np.random.random(size=(Natoms, 3))*10
	#zero velocity
	velocity=[[np.random.random(),  np.random.random(),  np.random.random()] for i in range(Natoms)]
	#print(get_acceleration(coord,mass))
	#coord, velocity=verlet(coord,  mass,  velocity,Nsteps=100,dt=0.1)
	coord1,velocity1=leap_frog(coord,  mass,  velocity,Nsteps=10,dt=0.05)
	os.system("move movie1.xyz movie_forward.xyz")
	leap_frog(coord1,mass,-1*velocity1,Nsteps=10,dt=0.05)
	os.system("move movie1.xyz movie_backward.xyz")

if __name__ == '__main__':
	main()