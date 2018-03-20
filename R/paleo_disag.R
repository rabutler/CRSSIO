
#this is a direct resampling of LB tribs - no scale factor applied
rm(list=ls())

nsim=1 # number of simulations	


l=1244# length of each simulation (yrs)

nsite=29

#dataD= matrix(scan("input/LF_total.txt"), ncol=20, byrow=TRUE) # obs data to be used for disag (total flow) 20 cols, monthly data #(20x1212)

data  =matrix(scan("input/CRB_CY_MON_IV_29_06_08.txt"), ncol=29, byrow=T) #intervening natural flow data - monthly CY text file


an_flw= matrix(scan("input/LF_06_08.txt"), ncol=2, byrow=TRUE) # observed annual flow for picking analog disag yr

al=length(an_flw)


len=length(data[,1])/12 # how many yrs of observed data


dat_a= array(data=NA,dim=c(len,12,nsite)) # matrix for observed values - row(yr), col (month), index (site)

pp=4



disag=array(data=NA,dim=c(l,12,nsite,nsim)) #matrix for disag values - row(yr), col (month), index1 (site), index2(sim #)

index_mat=matrix(ncol=nsim, nrow=l) # matrix for recording yr index for disag (optional)

sf_mat= matrix(ncol=nsim, nrow=l)#matrix for recording scale factors used (optional)

#this loop moves observed monthly data from 2d matrix to 3d array

mgn=length(dat_a[1,1,])

for(j in 1:mgn){

s=1
e=12

	for(i in 1:len){

	dat_a[i,,j]=data[s:e,j]

	s=s+12
	e=e+12
	}

}

seqs= matrix(scan("input/Meko.txt"), ncol=2, byrow=T) #read in annual synthetic data for disag


temp=array(data=NA,dim=c(l,12,nsim)) # temporary matrix for storing time only disag


for(j in 1:nsim){

#this picks the 1st year for disag based only on the annual flow
	k=sqrt(len) #number of neighbors
	
	Flow=seqs[1,(j+1)]
	
	D=abs(an_flw[,2]-Flow)
	
	Delta=cbind(an_flw[,1],D) #combines difference and corresponding year into one matrix
	
	Delta_sort=cbind(Delta[,1][order(Delta[,2])], sort(Delta[,2])) #reorders the delta matrix based on distances
	
	kmatrix = Delta_sort[1:k,1:2] #selects the "k-nearest-neighbors" from Delta_sort 
	
	weight = matrix(nrow=k, ncol=1) # defines matrix for weights
 		
	rnk=rank(kmatrix[,2]) #ranks distances for purpose of generating weights
		
		for(i in 1:k){
	
			weight[i,1] = 1/(rnk[i]) #fills weighting matrix
	
		}

	z = sum(weight) # sums weights 

	weights = weight/z	#divides weights by sum of weights so cumulative probability = 1
		
	N=sample(kmatrix[,1], 1, replace = TRUE, prob=weights) #Selects a year to be "nearest neighbor"
		

	pos=N-(an_flw[1,1]-1) # index for selected yr

	SF=Flow/(an_flw[pos,2]) # scaling factor to apply for disag

	temp[1, ,j]=dat_a[pos,,mgn]*SF	

	index_mat[1,j]=N
	sf_mat[1,j]=SF

	
        disag[1, ,1:20,j]=dat_a[pos, ,1:20]*SF
		disag[1, ,21:29,j]=dat_a[pos, ,21:29]

# now that one year has been disaggregated, the remaining years in the trace use annual flow and also december of last yr (CY)

		for(h in 2:l){


		k=sqrt(len) #number of neighbors

		Flow=seqs[h,(j+1)]
		D=2:len
	
        		
			D=abs(an_flw[,2]-Flow) # annual as the only selection criteria
			Delta=cbind(an_flw[,1],D) #combines difference and corresponding year into one matrix #these use just m.a.f


			Delta_sort=cbind(Delta[,1][order(Delta[,2])], sort(Delta[,2])) #reorders the delta matrix based on distances
		
			kmatrix = Delta_sort[1:k,1:2] #selects the "k-nearest-neighbors" from Delta_sort 
	
			weight = matrix(nrow=k, ncol=1) # defines matrix for weights
 		
			rnk=rank(kmatrix[,2]) #ranks distances for purpose of generating weights
		
			for(i in 1:k){
	
			weight[i,1] = 1/(rnk[i]) #fills weighting matrix
	
			}

			z = sum(weight) # sums weights 
	
			weights = weight/z	#divides weights by sum of weights so cumulative probability = 1
		
			N=sample(kmatrix[,1], 1, replace = TRUE, prob=weights) #Selects a year to be "nearest neighbor"
		
			pos=N-(an_flw[1,1]-1) # index for selected yr
	
			SF=Flow/(an_flw[pos,2]) # scaling factor to apply for disag

			
       

			#disag[h, , ,j]=dat_a[pos, ,]*SF
       		disag[h, ,1:20,j]=dat_a[pos, ,1:20]*SF
		    disag[h, ,21:29,j]=dat_a[pos, ,21:29]
			}
			}
			


## output to "flat" file


disagmat=matrix(ncol=nsite, nrow=(l*12*nsim))

for(i in 1:nsite){
p=1

for(k in 1:nsim){


for(j in 1:l){
disagmat[p:(p+11),i]=disag[j,,i,k]
p=p+12
}
}
}

write(t(disagmat), file="output/paleo_disag_4crss.txt", ncol=nsite)

#seqs=0

#source("code\\plottingCode\\plot_all.rm")




