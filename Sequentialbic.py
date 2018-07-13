import csv
import numpy as np
import pandas as pd

count = 0

# dataset 

f = open('/home/priyankit/Desktop/PROJ/Biclustering of Expression Data with evolutionary computing/yeast_dataset(5x5).csv','r')
csvR = pd.read_csv(f) 		# DataFrame
wk = csvR.values     		# numpy array and contains float values in rows AND rows(genes)-29 columns(conditions)-17 
x,y = wk.shape	     		# Dimensions of data

# genetic operators

no_of_gen = 10			# Number of Generations
no_of_bic = 5			# Number of Biclusters
wv = 1		
wrow = 1			# Weight of rows
wcol = 10			# Weight of columns
popsize = 4    		# Population Size
thresold = 1000		
pc = 0.85			# crossover probability
pm = 0.2			# mutation probability
elitism = 0.90			# elitism

index = 0 
population = np.array([])	# population - numpy array

# Defining blank numpy arrays

fitnessvalue = np.array([])	# fitness value - numpy array
msr = np.array([])		# mean squared residue - numpy array
var = np.array([])		# variance - numpy array
no_of_rows = np.array([])
no_of_col = np.array([])
volume = np.array([])		# volume - numpy array

avg_msr = np.array([])		# average mean squared residue - numpy array
avg_var = np.array([])		# average variance - numpy array
avg_vol = np.array([])		# average volume - numpy array
avg_fitness = np.array([])	# average fitness - numpy array

best_msr = np.array([])		# best mean squared residue - numpy array
best_var = np.array([])		# best variance - numpy array
best_vol = np.array([])		# best volume - numpy array
best_fitness = np.array([])	# best fitness - numpy array


final_msr = np.array([])	# final mean squared residue - numpy array
final_fitness = np.array([])	# final fitness - numpy array
final_variance = np.array([])	# final variance - numpy array
final_rows = np.array([])
final_col = np.array([])
final_vol = np.array([])	# final volume - numpy array


#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- N O R M A L I Z A T I O N   F U N C T I O N ------------------------------------
#----------------------------------------------------------------------------------------------------------------------

def normalize(x):		# x is a numpy array
	return ((np.array(x)- np.min(x))/(np.max(x) - np.min(x)))

#----------------------------------------------------------------------------------------------------------------------	
#--------------------------------------------E X P R E S S I O N   M A T R I X ----------------------------------------
#----------------------------------------------------------------------------------------------------------------------

# < N*O*T*E > : wk is the initial expression matrix
# < N*O*T*E > : expression_matrix is the normalized expression matrix

expression_matrix = np.array(map(normalize,wk))

#cnt = 0				# for counting through rows of expression_matrix
#for i in wk:
	#print i		# i is of type numpy array
#	expression_matrix = np.insert(expression_matrix,cnt,normalize(i))
#	cnt = cnt + 1

expression_matrix = np.array(expression_matrix).reshape(x,y) 	# As normalization is returning single vector only

row_total , col_total = expression_matrix.shape
stringlength = row_total + col_total

parent = np.array([0]*(popsize*stringlength)).reshape(popsize,stringlength)
offspring = np.array([0]*(popsize*stringlength)).reshape(popsize,stringlength)
results = np.array([0]*(no_of_bic*stringlength)).reshape(no_of_bic,stringlength)
em_weight = np.array([0]*(row_total*col_total)).reshape(row_total,col_total)


#-----------------------------------------------------------------------------------------------------------------------
#--------------------------------------- F I T N E S S  F U N C T I O N ------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

def fitness(x):
	global count
	global index
	global msr
	global var
	global volume
	global fitnessvalue
	global no_of_rows
	global no_of_col
	
	count1 = 0
  	count2 = 0
	penalty = 0
  	A = np.array([])
  	B = np.array([])
	x = np.array(x.tolist())

	# Finding a unique set of genes and conditions 
	for j in range(0,stringlength):
		if (x[j]==1) and (j < row_total):
  			A = np.insert(A,count1,j)
  			count1 = count1 + 1
  		if (x[j]==1) and (j >= row_total):
			B = np.insert(B,count2,j - row_total)
	  		count2 = count2 + 1

	#print "A"
	#print A
	#print "B"
	#print B

	submatrix = np.array([0.0]*(A.shape[0] * B.shape[0])).reshape(A.shape[0],B.shape[0]) # A Bicluster corresponds to submatrix that exhibits some coherant tendency
	
	#	BICLUSTERING
	if A.shape[0]==0 or B.shape[0]==0 :
		
		#print "index "+str(index)
		variance = 0
		msr = np.insert(msr,index,0.0)
		#print "\nmsr\n"
		#print msr
		residue = 0
		var = np.insert(var,index,0.0)
		#print "\nvar\n"
		#print var
		w_d = 0
		penalty = 0
		#print "0x0"
		volume = np.insert(volume,index,0.0)
		#print "\nvolume\n"
		#print volume
		fitnessvalue =	 np.insert(fitnessvalue,index,800)
		#print "\nfitnessvalue\n"
		#print fitnessvalue
		no_of_rows = np.insert(no_of_rows,index,0)
		#print "\nno of rows\n"
		#print no_of_rows
		no_of_col = np.insert(no_of_col,index,0)
		#print "\nno of col\n"
		#print no_of_col
		index = index + 1
		
		
	else:
		for i in range(0,A.shape[0]):
	  		for j in range(0,B.shape[0]):
				submatrix[i,j] = expression_matrix[int(A[i]),int(B[j])]
				penalty = penalty + em_weight[int(A[i]),int(B[j])]
	
		submat = pd.DataFrame(submatrix)			# BICLUSTER
		#print "\n---SUBMAT---\n"
		#print submat
	
								# DEFINATION <1>
		row_mean = np.array(submat.apply(np.mean,axis=1))	# ROW MEAN OF BICLUSTER or DEFINATION OF BASE OF THE GENE
		#print "\nROW MEAN\n"
		#print row_mean  	
		col_mean = np.array(submat.apply(np.mean,axis=0))	# COLUMN MEAN OF BICLUSTER or DEFINATION OF BASE OF THE CONDITIONS
		#print "\nCOLUMN MEAN\n"
		#print col_mean  	
		mat_mean = np.array(pd.DataFrame(np.mean(submat)).apply(np.mean,axis=0)) # BASE OF THE BICLUSTER i.e. THE MEAN OF ALL ENTRIES IN BICLUSTER
		#print "\nMATRIX MEAN\n"
		#print mat_mean	
		newmat_residue = np.array([0.0]* submat.shape[0] * submat.shape[1]).reshape(submat.shape[0],submat.shape[1])
		newmat_msr =     np.array([0.0]* submat.shape[0] * submat.shape[1]).reshape(submat.shape[0],submat.shape[1])	
	
		# RESIDUE OF AN ENTRY e(i,j) OF A BICLUSTER(I,J) OR DEFINATION <2>
		for i in range(0,newmat_residue.shape[0]) :
 			for j in range(0,newmat_residue.shape[1]):
  				newmat_residue[i,j] = ( submatrix[i,j] - row_mean[i] - col_mean[j] +mat_mean )  # submatrix is same as submat
				newmat_msr[i,j] = newmat_residue[i,j] ** 2
		#print "\n---NEWMAT_RESIDUE---\n"
		#print pd.DataFrame(newmat_residue)
	
		#print "\n---NEWMAT_MSR---\n"
		#print pd.DataFrame(newmat_msr)
	
		# < N*O*T*E > : RESIDUE IS AN INDICATOR OF THE DEGREE OF COHERANCE OF AN ELEMENT W.R.T REMAINING ONES IN THE BICLUSTER
		# < N*O*T*E > : LOWER THE RESIDUE STRONGER THE COHERANCE
		residue = (np.mean(newmat_residue))*10**4								# PREV : (np.mean(newmat_residue))*10**4
		#print "\nresidue : "+str(residue)
	
		msr = np.insert(msr,index,(np.mean(newmat_msr))/(newmat_msr.shape[0] * newmat_msr.shape[1]))		# DEFINATION <3> 
		#print "\n\n\nM S R	MATRIX\n"	
		#print msr	
		row_var = np.array(submat.apply(np.var,axis=1))
		#print "\n\nrow_var\n\n"
		#print row_var
		variance = (sum(row_var))*10**3
		#print "\n\nVariance\n\n"
		#print variance
		# < N*O*T*E > : USING VARIANCE WE WANT TO GUARENTEE THAT THE BICLUSTER CAPTURESGENES EXHIBITING FLUCTUATING YET COHERANT TRENDS UNDER SOME 				SET CONDITIONS 
		var = np.insert(var,index,variance)
		#print "\n\nVAR MATRIX\n\n"
		#print var

		w_d = wv * ( wrow * (thresold/A.shape[0]) + wcol * ( thresold / B.shape[0] )  )
		#print "\nw_d : "+str(w_d)		
		
		fitnessvalue = np.insert(fitnessvalue,index, (residue/thresold)+(1/variance)+w_d+penalty )
					
		#print "\n\nFITNESS VALUE\n\n"
		#print fitnessvalue
		#print count
		#count = count + 1
		volume = np.insert(volume,index,w_d)
		#print "\n\nVOLUME\n\n"
		#print volume
	  	no_of_rows = np.insert(no_of_rows,index,count1 - 1)
		#print "\n\nno_of_rows\n\n"
		#print no_of_rows
		no_of_col = np.insert(no_of_col,index,count2 - 1)
		#print "\n\nno of col\n\n"
		#print no_of_col		
		index = index + 1
		#print "--------------------------------------------------------------------------"
		
	
	
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------- P A R E N T  S E L E C T I O N  F U N C T I O N ----------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

# < N*O*T*E > : Parent selection is bsed on tournament selection as per paper . To check whether the folowing is as per paper
def selection_parents(population) :
	global parent
	count = 0
	while count < popsize :
		i = np.random.choice(range(0,2*popsize),size=1)[0]
		j = np.random.choice(range(0,2*popsize),size=1)[0]
		#i = np.random.choice(range(0,popsize),size=1)[0]
		#j = np.random.choice(range(0,popsize),size=1)[0]
		#print i
		#print j
		if i!= j :
			if np.isnan(fitnessvalue[j]) == False  and np.isnan(fitnessvalue[i]) == False :
				if fitnessvalue[i] < fitnessvalue[j] :
					if np.random.uniform(0,1,1)[0] < elitism :
						parent[count] = population[i]
				else:
					if np.random.uniform(0,1,1)[0] < elitism :
						parent[count] = population[j]
		else:
			parent[count] = population[i]
		count = count + 1

	#print "\n\nPARENT\n\n"
	#print parent


#-------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------- C R O S S O V E R  F U N C T I O N --------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

# < N*O*T*E > : crossover swaps some genetic material between two or more individuals
def crossover() :
	global offspring
	global population
	global parent 
	s1 = np.array([])
	s2 = np.array([])
	part = np.array([])

	for i in list(range(0,popsize-1,2)) :
		crosspoint = np.random.choice(range(0,stringlength),size=1)[0]
		#print "\ni\n"
		#print i
		#print "crosspoint"
		#print crosspoint
		s1 = parent[i]
		#print "\ns1\n"
		 	#print s1
		s2 = parent[i+1]
		#print "\ns2\n"
		#print s2
		#print "\ns1[crosspoint:stringlength]\n"
		#print s1[crosspoint:stringlength]
		#print "\ns2[crosspoint:stringlength]\n"
		#print s2[crosspoint:stringlength]
		part = s1[crosspoint:stringlength]
		s1[crosspoint:stringlength] = s2[crosspoint:stringlength]
		s2[crosspoint:stringlength] = part
		#print "\ns1\n"
		#print s1
		#print "\ns2\n"
		#print s2
		offspring[i] = s1
		offspring[i+1] = s2
	
	#print "\n\nCrossover popsize\n"	
	#print parent[0:popsize]
	population[0:popsize] = parent[0:popsize]
	population[popsize:(2*popsize)] = offspring[0:popsize]



#-------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------  M U T A T I O N   F U N C T I O N  ---------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

# < N*O*T*E > : mutation changes a small part of the genetic material of an individual to a new random value
def mutation() :
	global population
	for i in list(range(0,popsize-1,2)) :
		mut = int(2 * popsize * stringlength * 0.01)
		for iter in range(0,mut) :
			c = np.random.choice(range(0,stringlength),size=1)[0]
			r = np.random.	choice(range(0,popsize),size=1)[0]
			if population[r,c]==1:
				population[r,c] = 0
			else:
				population[r,c] = 1
				
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------   E B I  F U N C T I O N   --------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

def ebi(expression_matrix,thresold):
	global index
	global population
	global avg_fitness
	global avg_msr
	global avg_var
	global avg_vol
	global best_fitness
	global best_msr
	global best_var
	global best_vol
	global em_weight

	max_iter = 0
	population = np.array([0]*(popsize*2*stringlength)).reshape(popsize*2,stringlength) # 40 row n 46 col

	for i in range(0,popsize):	
		for j in range(0,stringlength):
			population[i,j] = np.random.choice([0,1],size=1)[0]
	#population = np.array(np.matrix([[0,0,1,0,0,0,0,0,0,1],[0,0,1,1,0,1,0,1,0,1],[1,1,0,0,1,0,0,1,0,1],[1,0,0,0,1,1,1,0,1,1],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0]]))
	#print "\n\nPOPULATION\n\n" 		# CORRECT WAY TO SHOW POPULATION
	#for i in population:
	#	print i
	#	print "\n\n"
	
	index = 0
	map(fitness, population)
	
	
	#popDF = pd.DataFrame(population)
	#tPopulation = popDF.apply(fitness,axis=1)
	#print "tPopulation"
	#print popDF.apply(fitness,axis=1)
	
	while max_iter < no_of_gen :
		#print "=================================================================================="
		#print "=================================================================================="
		selection_parents(population)	# 20 rows and 46 columns
		#print "\nPARENT\n"
		#print parent
		
			
		if np.random.uniform(0,1,1)[0]  < pc :
			crossover()
			#print "\n*************** AFTER CROSSOVER ****************\n"
			#print "\nPOPULATION\n"
			#print population
			#print "\nOFFSPRING\n"
			#print offspring
		if np.random.uniform(0,1,1)[0]  < pm :
			mutation()
			#print "\n*************** AFTER MUTATION ****************\n"
			#print "\nPOPULATION\n"
			#print population
		
		index = 0
		

		map(fitness, population)
			
		#tPopulation = popDF.apply(fitness,axis=1)
		
		#print "\n\nGen "+str(max_iter)
		avg_fitness = np.insert(avg_fitness,max_iter,fitnessvalue.mean())
		#print "\n\navg_fitness\n\n"
		#print fitnessvalue.mean()
		avg_msr = np.insert(avg_msr,max_iter,msr.mean())
		#print "\n\navg_msr\n\n"	
		#print msr.mean()
		avg_var = np.insert(avg_var,max_iter,var.mean()) 
		avg_vol = np.insert(avg_vol,max_iter,volume.mean())
		#max_iter = max_iter+1
		
		best_fitness = np.insert(best_fitness,max_iter,np.array(fitnessvalue).min())
		best_msr = np.insert(best_msr,max_iter,np.array(msr).min())
		best_var = np.insert(best_var,max_iter,np.array(var).max())
		best_vol = np.insert(best_vol,max_iter,np.array(volume).max())
		max_iter = max_iter+1
		
		
	index_min = np.nanargmin(fitnessvalue)
	#print "\n\n\n\n----------------------------------------\nindex_min"
	#print index_min
	#print fitnessvalue.shape
	#print population.shape
	best_ind = population[index_min]
	
	if msr[index_min] < thresold - 1 :
		mylist = {"index" : index_min , "best_ind" : best_ind}		# DICTIONARY
		#print "\n\n--------------------------------- msr[index_min] -------------------------------------------------------\n\n"
		#print msr
		
		count1 = 0
		count2 = 0
		A = np.array([])
		B = np.array([])
		for j in range(0,best_ind.shape[0]) :
			if ( best_ind[j] == 1 ) and ( j < row_total) :
				A = np.insert(A,count1,j)
	    			count1 = count1 - 1
			if ( best_ind[j] == 1 ) and ( j >=row_total ) :	
				B = np.insert(B,count2,j-row_total)
	    			count2 = count2 - 1
		#print "\n\n A - EBI FUNCTION---------------------------------------------------------------------------------------\n\n"
		#print A
		#print "\n\n B - EBI FUNCTION\n\n"
		#print B
		for i in range(0,A.shape[0]) :
			  for j in range(0,B.shape[0]) :
	    			em_weight[A[i],B[j]] = em_weight[A[i] , B[j]] + 1
		return mylist
	
	
#---------------------------------------------------------------------------------------------------------------------------------
#------------------------------  S E Q U E N T I A L   C O V E R I N G   F U N C T I O N -----------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------

def sequential_covering(expression_matrix):
	global results
	global final_fitness
	global final_msr
	global final_variance
	global final_rows
	global final_col
	global final_vol

	fr = np.array([0]*(no_of_bic*6)).reshape(no_of_bic,6) # 6 columns for (1)fitness (2)row variance (3)msr (4)volume (5)no of rows (6)no of columns
	cols = ["--FITNESS-- ","ROW VARIANCE","---MSR---","VOLUME","NO OF ROWS","NO OF COLUMNS"]
	names = [str(unichr(65+i)) for i in range(0,no_of_bic)] # aplphabetical index for rows or entities
	final_results = pd.DataFrame(fr,index=names,columns=cols)
	#print final_results

	res_count = 0

   	iter = 0
	
	#ebi(expression_matrix,thresold) # testing purpose to be contd from here on & PREV thresold was 1000 hardcoded	{ TO BE CHECKED }

    	while iter < no_of_bic :	
   		rlist = ebi(expression_matrix,thresold)
		#print "\n\n==================================================================================\n\n"
		#print rlist
		#print "\n\n==================================================================================\n\n"
      		if len(rlist)!=0:
	      		results[res_count] = rlist["best_ind"]
	      		index_var_msr = rlist["index"]
			
	      		final_fitness = np.insert(final_fitness,res_count,fitnessvalue[index_var_msr])
			#print final_fitness
	      		final_msr = np.insert(final_msr,res_count,msr[index_var_msr])
			#print final_msr
	      		final_variance = np.insert(final_variance,res_count,var[index_var_msr])
	      		final_rows = np.insert(final_rows,res_count,no_of_rows[index_var_msr])
	      		final_col = np.insert(final_col,res_count,no_of_col[index_var_msr])
	      		final_vol = np.insert(final_vol,res_count,volume[index_var_msr])
	      	    
	      		res_count = res_count+1
	      		iter = iter+1
   		else:
    			continue
	
	#print results
	
  	final_results = pd.concat([ pd.DataFrame(final_fitness[0:no_of_bic]) , pd.DataFrame(final_variance[0:no_of_bic]) , pd.DataFrame(final_msr[0:no_of_bic]) , pd.DataFrame(final_vol[0:no_of_bic]) , pd.DataFrame(final_col[0:no_of_bic]) , pd.DataFrame(final_rows[0:no_of_bic]) ] , axis=1)
  
        print final_results
	print avg_fitness[0:no_of_gen] 
 	print avg_vol[0:no_of_gen]
	print avg_var[0:no_of_gen]
	print avg_msr[0:no_of_gen]


#--------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	sequential_covering(expression_matrix)



# TO CHECK
# algo of parent selection
# algo of crossover
# algo of mutation

