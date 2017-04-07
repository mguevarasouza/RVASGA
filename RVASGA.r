
##############################################################
#RVASGA created by Mauricio Guevara March 2015
##############################################################

#Function to run the GA by chromosome or in the entre genome
#############################################################
RVASGA<-function(win, chromosome){
  library(WriteXLS)
  library("hash")
  if(missing(chromosome)){
    runAllGA(win)
  }
  else{
    runChromoGA(win, chromosome)
  }

}
#############################################################


#############################################################
#Main function to run the GA in one choromosome and store results
#in excel Param win is used as an offset to determine the size of
#the region where genes with effect are going to be searched
##############################################################
runChromoGA<-function(win, chromosome){
  print("GA by Chromosome")
  dt = matrix()
  cv = matrix()
  off<<-win
  cached<<-0
  testDir<-getwd()
  start.time <- Sys.time()
  print(paste0("Starting analyses of chromosome: ", chromosome))
  hashTable<<-hash()
  dirToSearch<-paste(testDir, sep="/",chromosome)
  setwd(dirToSearch)
  tryGA()
  s<-getFitness(lg2)
  a<- c(chromosome,paste(lg2, collapse = " "),s)
  dt <- rbind(dt,paste(a,collapse= ","))
  df <- as.data.frame(dt)
  setwd(testDir)
  # Write CSV in R
  name<-paste0(chromosome,".csv")
  write.table(df, file=name, sep="",col.names = F, row.names = F, na ="", quote = F)

  dt2 = matrix()
  for(i in lg2){
   ft <-getFitness(i)
   a2<-c(i, ft)
   dt2 <- rbind(dt2,paste(a2,collapse= ","))
  }

  df2<- as.data.frame(dt2)
  write.table(df2, file=name, append = TRUE, sep="  ",col.names = F, row.names = F, na ="", quote = F)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

#############################################################
#Main function to run the GA in all chromosomes and store results in excel
#Param win is used as an offset to determine the size of the
#region where genes with effect are going to be searched
##############################################################
runAllGA<-function(win){
  print("Whole Genome analyses")
  dt = matrix()
  cv = matrix()
  off<<-win
  cached<<-0
  testDir<-getwd()
  start.time <- Sys.time()
  for(j in 1: 22) {
    print(paste0("Starting analyses of chromosome: ", j))
    hashTable<<-hash()
    dirToSearch<-paste(testDir, sep="/",j)
    setwd(dirToSearch)
    tryGA()
    s<-getFitness(lg2)
    a<- c(j,lg2,s)
   dt <- rbind(dt,paste(a,collapse= ","))
  }
  df <- as.data.frame(dt)
  setwd(testDir)
  #Write to Excel
  WriteXLS("df", ExcelFileName = "Genome.xls", SheetNames = NULL, perl = "perl")
  # Write CSV in R
  write.table(df, file="Genome.csv", sep=",   ")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  print(cached)
}


##############################################################
#Function that reads the genes with effect
##############################################################
readGenesWithEffect<-function(fileEffectGenes){
  effect<-matrix(nrow=0,ncol=2)
  if(file.exists(fileEffectGenes)){
    mydata=read.csv(fileEffectGenes, header=FALSE)
    for(i in 1:length(mydata[,1])){
	  effect<-rbind(effect,mydata[i,])
	}
  }
effect
}

##############################################################
#Functions to prepare the genetic algorithm
##############################################################

##############################################################
#Function that gets the starting position of the gene
##############################################################
getStart<-function(fileName){
	con<-file(fileName, "r")
	line<-read.table(con,nrow=1)
	id<-line[1]
	str<-paste(unlist(id), collapse="")
 	start<-unlist(strsplit(str,":"))[3]
	close(con)
	start
}

##############################################################
#Function that calculates the chromosome length
##############################################################
calcCL<-function(){
	first<<-genesMat[1,2]
	last<<-genesMat[dim(genesMat)[1],2]
	cl<<-ceiling(log2(last-first))
}

##############################################################
#Function that reads all the files in a folder and makes a table
# with the starting point of each gene and loads in memory all
#the genes
##############################################################
readGenes<-function(){
  library(SKAT)
  library(GA)
  options(warn=-1)
  print("Constructing gene matrix")
  filesMat<-list.files(path = ".", pattern = "_geno.txt", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  t<-c()
  for(i in 1:length(filesMat)){
    assign(gsub("_geno.txt", "", filesMat[i]),read.table(filesMat[i]) ,envir=.GlobalEnv)
    t<-c(t,getStart(filesMat[i]))
  }

  calcObj(filesMat)
  genesMat<<-data.frame(filesMat)
  genesMat<<-cbind(genesMat,as.numeric(t))
  genesMat<<-genesMat[order(genesMat[,2]),]
  remove(filesMat)
  remove(t)
  calcCL()
}

##############################################################
#function that return the pvalue of a list of genes
##############################################################
getFitness<-function(listGenes){
  if(length(listGenes) > 0){
    assign("file", listGenes[1])

	if(length(listGenes) > 1){
	  tmp<-"files"
		for(i in 2:length(listGenes)){
		  tmp<-paste(tmp,i,sep="")
		  assign(tmp, listGenes[i])
		  tmp<-"files"
		}
	}
	a<-get(file)
	a$V1<-NULL
	comb<-(t(a))
	files<-ls(pattern="files*")
	for (f in files) {
	  if (!is.null(f)) {
		b<-get(get(f))
		b$V1<-NULL
		comb<-cbind(comb,t(b))
	  }
	}
	s<-SKAT(comb, objT)$p.value

	}else{
	  1
	}
}

##############################################################
#Function that calculates the fitness of a given chromosome
#In this case chromosome refers to a GA chromosome
##############################################################
fitness<-function(chromosome){
  rng<-calcRange(chromosome,off)
  listGenes<-getGN(rng)
  if(length(listGenes) > 0){
    listGenesk<-paste(listGenes,collapse=" ")
    if(has.key(listGenesk,hashTable)){
	  s<-hashTable[[listGenesk]]
	  cached<<-cached+1
	}else{
	  s<-getFitness(listGenes)
	  .set(hashTable, listGenesk,s)
	}
  -s
  }else{
    -1
  }
}

##############################################################
#Fitness that calculates the range of the search
##############################################################
calcRange<-function(chromosome, off){
  a<-paste(chromosome, collapse='')
  b<-strtoi(a, base = 2)
  g <- round(b * last / 2^cl )
  LB<- (g-off)
  UB<- (g+off)
  if(LB < first){
	LB<-first
  }
  if(UB > last){
	UB<-last
  }
  r<-cbind(LB, UB)
}

##############################################################
#Function that prepares the SKAT model
##############################################################
calcObj<-function(filesMat) {
   a<-gsub("_geno.txt", "", filesMat[1])
   p<-read.table(paste(a,"_pheno.txt",sep=""),header=F)
  objT<<-SKAT_Null_Model(p$V2 ~ 1, out_type="D")
}

##############################################################
#Function that gets the list of genes depending on the range
##############################################################
getGN<-function(rng){
  LB<-rng[1]
  UB<-rng[2]
  LI<-0
  UI<-0
  A<-TRUE
  B<-TRUE
  for(i in 1:length(genesMat[,1])){
	if(genesMat[i,2] >= LB && A){
      LI<-i
	  if(LI > 1){
		LI<-LI-1
	  }
	  A<-FALSE
	}
	if(genesMat[i,2] >= UB && B){
	  UI<-(i-1)
	  B<-FALSE
	  break
	}
  }
  H<-genesMat[LI:UI,1]
  H<-levels(H)[H]
  H<-gsub("_geno.txt", "", H)
}

############################################################
#Functions to create folders and move files
############################################################
#Move file to the correct folder
moveToFolder<-function(fileName) {
  con<-file(fileName, "r")
  line<-read.table(con,nrow=1)
  id<-line[1]
  str<-paste(unlist(id), collapse="")
  chromosome<-unlist(strsplit(str,":"))[2]
  from<-fileName
  to<- paste(chromosome,"/",fileName, sep='')
  close(con)
  my.file.rename(from,to)
  from<-gsub("geno", "pheno", from)
  to<-gsub("geno", "pheno", to)
  my.file.rename(from,to)
  from<-gsub("pheno", "mapping", from)
  to<-gsub("pheno", "mapping", to)
  my.file.rename(from,to)
}

##############################################################
#Function to create a folder for each chromosome and move the
#files to the corresponding folder
##############################################################
splitFiles<-function(){
  filesMat<-list.files(path = ".", pattern = "_geno.txt", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

  for(i in 1:length(filesMat)){
	 moveToFolder(filesMat[i])
	}
}

##############################################################
#Utility function to move the file
##############################################################
my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
	file.rename(from = from,  to = to)
}

##############################################################
#Functions to get the results after the GA run
##############################################################

##############################################################
#function to print genes of a population
##############################################################
printGenes<-function(p){
d<-NULL
 for(i in 1:length(p[,1])){
  rng<-calcRange(p[i,],0)
  listGenes<-getGN(rng)
  d<-rbind(d,listGenes)
 }
 print(unique(d))
}

##############################################################
#function that gets the best solution of the GA
##############################################################
getBestGenes<-function(kb){
  pivot<-unlist(lapply(kb[length(kb)],tail,1))
  rng<-calcRange(pivot,off)
  listGenes<-getGN(rng)
  listGenes
}

##############################################################
#function that gets the best solution of the second GA
##############################################################
getBestGenes2<-function(kb){
  ch<-unlist(lapply(kb[length(kb)],tail,1))
  genes<-getSelectedGenes(ch)
}

################################################
#Functions to run the GA
#################################################

##############################################################
#Function to run GA with one offset and keeps the best solution
##############################################################
tryGA<-function(){
  readGenes()
  clear(hashTable)
  print("Starting Genetic Algorithm 1");
  GA<-ga(type="binary", fitness = fitness, pcrossover=.6, pmutation=.2, elitism =0, nBits = cl, run = 2, maxiter=50, popSize=200, keepBest = TRUE)
  lg<<-getBestGenes(GA@bestSol)
  calcCL2(lg)
  clear(hashTable)
  print("Starting Genetic Algorithm 2");
  GA2<-ga(type="binary", fitness = fitness2, pcrossover=.6, pmutation=.2, elitism =0, nBits = cl2, run = 5, maxiter=100, popSize=100, keepBest = TRUE)
  lg2<<-getBestGenes2(GA2@bestSol)
}

##############################################################
#functions for the second genetic algorithm
#lg is the list of genes yielded by the first GA
##############################################################

##############################################################
#Function to calculate the chromosome length
##############################################################
calcCL2<-function(listGenes){
  cl2<<-length(listGenes)
}


##############################################################
#Get the genes with index 1 in chromosome
##############################################################
getSelectedGenes<-function(ch){
  idx<-which(ch==1)
  t<-lg[idx]
  t
}

##############################################################
##fitness function for the second genetic algorithm,
#Returns the fitness of a chromosome
##############################################################
fitness2<-function(ch){
  listGenes<-getSelectedGenes(ch)
  if(length(listGenes) > 0){
	listGenesk<-paste(listGenes,collapse=" ")
	if(has.key(listGenesk,hashTable)){
	  s<-hashTable[[listGenesk]]
	  cached<<-cached+1
	}else{
	  s<-getFitness(listGenes)
	  .set(hashTable, listGenesk,s)
	}
	-s
  }else{
	-1
  }
}
