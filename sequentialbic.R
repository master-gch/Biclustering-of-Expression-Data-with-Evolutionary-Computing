    library(XLConnect)
    library(matrixStats) 
    library(graphics)
    
    wk <- read.csv('/home/priyankit/Desktop/PROJ/Biclustering of Expression Data with evolutionary computing/yeast_dataset(5x5).csv') 
    #print(wk)
    count=0
    no_of_gen<- 2
    no_of_bic<- 5
    wv<-1
    wrow<-1
    wcol<-10
    popsize<-4
    threshold<-1000
    pc<-0.95
    pm<-0.20
    elitism<-0.90
    
    fitnessvalue<-c()
    msr<-c()
    var<-c()
    no_of_rows<-c()
    no_of_col<-c()
    volume<-c()
    
    avg_msr<-c()
    avg_var<-c()
    avg_vol<-c()
    avg_fitness<-c()
    
    best_msr<-c()
    best_var<-c()
    best_vol<-c()
    best_fitness<-c()
    
    final_msr<-c()
    final_fitness<-c()
    final_variance<-c()
    final_rows<-c()
    final_col<-c()
    final_vol<-c()
    
    
    normalize <- function(x) 
    {
     return ((x - min(x)) / (max(x) - min(x)))
    }
    
    expression_matrix <- as.matrix(apply(wk,1 ,normalize))
    #print(expression_matrix)
    row_total<-nrow(expression_matrix)
    col_total<-ncol(expression_matrix)
    stringlength<-row_total+col_total
    #print(paste0("total rows",row_total))
    #print(paste0("total columns",col_total))
    parent<-matrix(0,nrow=popsize,ncol=stringlength)
    offspring<-matrix(0,nrow=popsize,ncol=stringlength)
    results<-matrix(0,nrow=no_of_bic,ncol=stringlength)
    em_weight<-matrix(0,nrow=row_total,ncol=col_total)
    
    
    fitness<-function(x)
    {
      count1<-1
      count2<-1
      penalty<-0
      A<-c()
      B<-c()
      	for(j in 1:stringlength)
      	{
      		if((x[j]==1) && (j<=row_total))
      		{
      		A[count1]<-j
      		count1<-count1+1
      		}
      		if((x[j]==1) && (j>row_total))
      		{
      		B[count2]<-j-row_total
      		count2<-count2+1
      		}
      	}
      #print("A\n")
      #print(A)
      #print("\n\nB\n")
      #print(B)
      submatrix<-data.frame()
      if(length(A)==0 | length(B)==0)
      {submatrix<-data.frame()
          variance<-0
          msr[index]<<-0
          #print("msr")
          #print(msr)
          residue<-0
          var[index]<<-0
          #print("var")
          #print(var)
          w_d<-0
          penalty<-0
          #print(paste0("0x0"))
          volume[index]<<-0
          #print("volume")
          #print(volume)
          fitnessvalue[index]<<-800
          #print("fitnessvalue")
          #print(fitnessvalue)
          no_of_rows[index]<<-0
          #print("no of rows")
          #print(no_of_rows)
          no_of_col[index]<<-0
          #print("no of col")
          #print(no_of_col)
          index<<-index+1
      }
      else
      {
      	for(i in 1:length(A))
      	{
      		for(j in 1:length(B))
      		{
      		submatrix[i,j]<-expression_matrix[A[i],B[j]]
      		penalty<-penalty+em_weight[A[i],B[j]]
      		}
      	}
      submat<-as.matrix(submatrix)
       #print("submat")
       #print(submat)
      row_mean<-apply(submat,1,mean)
      col_mean<-apply(submat,2,mean)
      mat_mean<-mean(submat)
      newmat_residue<-matrix(0,nrow=nrow(submat),ncol=ncol(submat))
      newmat_msr <- matrix(0,nrow=nrow(submat),ncol=ncol(submat))
      	for(i in 1:nrow(newmat_residue))
      	{
      		for(j in 1:ncol(newmat_residue))
      		{
      		newmat_residue[i,j]<-(submat[i,j]-row_mean[i]-col_mean[j]+mat_mean)
      		newmat_msr[i,j]<- newmat_residue[i,j]**2
      		}
      	}
      #print("newmat_residue")
      print(newmat_msr)
      residue<-(mean(newmat_residue))
      #print(residue)
      #msr[index]<<-residue
      msr[index]<<- sum(newmat_msr) / (dim(newmat_msr)[1] * dim(newmat_msr)[2])
      row_var<-rowVars(submat, rows=NULL, cols=NULL, NA.rm=TRUE)
      variance<-(sum(row_var))*10^3
      var[index]<<-variance
      w_d<-wv*(wrow*(threshold/length(A))+wcol*(threshold/length(B)))
      fitnessvalue[index]<<-(residue/threshold)+(1/variance)+w_d+penalty
      volume[index]<<-w_d
      no_of_rows[index]<<-count1-1
      no_of_col[index]<<-count2-1
      
      #print(count)
      count <<- count +1
      #print("w_d")
      #print(w_d)
      #print("fitnessvalue")
      #print(fitnessvalue)
      index<<-index+1
      }
      print("---------------------------------------------------------------------------")
    }
    
    selection_parents<-function(population)
    {
    count<-1
    	while(count<=popsize)
    	{
    	i<-sample(1:(popsize), 1)
    	j<-sample(1:(popsize), 1)
    	#print(i)
    	#print(j)
    		if(i!=j)
    		{
    			if(!is.na(fitnessvalue[j])  && !is.na(fitnessvalue[i]))
    			{
    				if(fitnessvalue[i]<fitnessvalue[j] )
    				{
    					if(runif(1,0,1)<0.90)
    					{
    					parent[count,]<<-population[i,]
    					}
    				}
    				else
    				{
    					if(runif(1,0,1)<0.90)
    					{
    					parent[count,]<<-population[j,]
    					}
    				}
    			}
    		}
  
    		else
    			parent[count,]<<-population[i,]
    		count<-count+1
  
    	}
    }
  
    crossover<-function(population)
    {
      #print(population)
    s1<-list()
    s2<-list()
    part<-list()
    	for(i in seq(1, popsize-1, 2))
    		{
    		crosspoint<-sample(1:stringlength,1)
    		s1<-parent[i,]
    		s2<-parent[i+1,]
    		part<-s1[crosspoint:stringlength]
    		s1[crosspoint:stringlength]<-s2[crosspoint:stringlength]
    		s2[crosspoint:stringlength]<-part
    		offspring[i,]<<-s1
    		offspring[i+1,]<<-s2
    	}
    #print("\n\nCrossover popsize\n")	
    #print(parent[1:popsize,])
    population[1:popsize,]<<-parent[1:popsize,]
    population[(popsize+1):(2*popsize),]<<-offspring[1:popsize,]
    
    #print("******************CROSSOVER***************")
    #print("POPULATION")
    #print(population)
    #print("OFFSPRING")
    #print(offspring)
    }
  
    mutation<-function(population)
    {
    	for(i in seq(1, popsize-1, 2))
    	{
    	mut<-as.integer(2*popsize*stringlength*0.01)
    		for(iter in 1:mut)
    		{
    		c<-sample(1:(stringlength),1)
    		r<-sample(1:(popsize),1)
    		if(population[r,c]==1)
    			population[r,c]<<-0
    		else
    			population[r,c]<<-1
    		}
    	}
      #print("******************MUTATION***************")
      #print("POPULATION")
      #print(population)
    }
  
    
    ebi<- function(expression_matrix,threshold)
    {
    max_iter<-1
    population <<- matrix(0, 2*popsize, stringlength)   
    	for(i in 1:popsize)
    	{
      		for(j in 1:stringlength)
    		{
        		population[i, j] <<- sample(c(0,1),1)
      	}
    	}
    #print("POPULATION")
    
    index<<-1
    population <<- matrix(c(0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0),nrow=8,ncol=10)
    #print(population)
    apply(population, 1, fitness)
  
     		while(max_iter<=no_of_gen)
     		{
     		 #print ("==================================================================================")
     		 #print ("==================================================================================")  
     		selection_parents(population)
     		#print(parent)
     		if(runif(1,0,1)<pc)
     			{
     		  crossover(population)
     		  #print("*******************AFTER CROSSOVER************************")
     		  #print("POPULATION")
     		  #print(population)
     		  #print("OFFSPRING")
     		  #print( offspring)
     		}
     		if(runif(1,0,1)<pm)
     		{
     		  mutation(population)
     		  #print("*********************AFTER MUTATION*************************")
     		  #print("POPULATION")
     		  #print(population)
     		}
     		index<<-1
     		apply(population,1,fitness)
     		
     		avg_fitness[max_iter]<<-mean(fitnessvalue)
     		avg_msr[max_iter]<<-mean(msr)
     		avg_var[max_iter]<<-mean(var)
     		avg_vol[max_iter]<<-mean(volume)
     		
     		best_fitness[max_iter]<<-min(fitnessvalue)
     		best_msr[max_iter]<<-min(msr)
     		best_var[max_iter]<<-max(var)
     		best_vol[max_iter]<<-max(volume)
     		max_iter<-max_iter+1
     		}
    
    index_min<-which.min(fitnessvalue)
    #print("index_min")
    #print(index_min)
    best_ind<-population[index_min,]
    
     if(msr[index_min]<threshold)
    {
      mylist<-list("index"=index_min,"best_ind"=best_ind)
      count1<-1
      count2<-1
      A<-c()
      B<-c()
      for(j in 1:length(best_ind))
      {
        if((best_ind[j]==1) && (j<=row_total))
        {
          A[count1]<-j
          count1<-count1-1
        }
        if((best_ind[j]==1) && (j>row_total))
        {
          B[count2]<-j-row_total
          count2<-count2-1
        }
      }
      #print(A)
      #print(B)
      for(i in 1:length(A))
      {
        for(j in 1:length(B))
        {
          em_weight[A[i],B[j]]<<-em_weight[A[i],B[j]]+1
        }
      }
      return(mylist)
    }
    
    }
    
    
    
     sequential_covering<-function(expression_matrix)
     {
       final_results<-as.data.frame(matrix(0,nrow=no_of_bic,ncol=6))
       colnames(final_results)<-list("--FITNESS-- ","ROW VARIANCE","---MSR---","VOLUME","NO OF ROWS","NO OF COLUMNS")
       res_count<-1
       iter<-1
     #  ebi(expression_matrix,threshold)
     
     		while(iter<=no_of_bic)
      		{
     		rlist<-ebi(expression_matrix,threshold)
     		print(rlist)
        		if(!length(rlist)==0)
        		{
       		results[res_count,]<<-rlist$best_ind
       		index_var_msr<-rlist$index
       		final_fitness[res_count]<<-fitnessvalue[index_var_msr]
       		print(final_fitness)
       		final_msr[res_count]<<-msr[index_var_msr]
       		final_variance[res_count]<<-var[index_var_msr]
       		final_rows[res_count]<<-no_of_rows[index_var_msr]
       		final_col[res_count]<<-no_of_col[index_var_msr]
       		final_vol[res_count]<<-volume[index_var_msr]


        		res_count<-res_count+1
        		iter<-iter+1
        		}
      		else
      		next
      		}
     #   print(results)
          for(i in 1:no_of_bic)
          {
          final_results[i,]<-cbind(final_fitness[i],final_variance[i],final_msr[i],final_vol[i],final_col[i],final_rows[i])
          }
      print(final_results)
      print(avg_fitness)
      print(avg_vol)
      print(avg_var)
      print(avg_msr)
    
    }
    sequential_covering(expression_matrix)
