## PharmaHacks 2019
## Team Novinger Inc.
## Robin, Ryan, Guang, Chris, Frederic


##### Challenge 1 ####

## Briefly, the task is to recreate the SMILE of a molecule based only on a set of given features.
## We break this down into 3 subtasks:
##   A) Given the exact molecular weight, and NOcount, we can deduce the chemical formula
##   B) Given the chemical formula, and a list of moieties, we can deduce possible bond graphs (aka. adjacency matrices)
##   C) Given a adjacency matrix, we can generate the SMILE and the 3-D structure


##### Setup #####

# install the following libraries if needed
library(hier.part) #https://stackoverflow.com/questions/17292091/rbinary-matrix-for-all-possible-unique-results
library(e1071) #http://ugrad.stat.ubc.ca/R/library/e1071/html/countpattern.html

library(partitions) #https://stackoverflow.com/questions/22218640/getting-all-combinations-which-sum-up-to-100-using-r
library(iterpc) #same link as above

library(readxl)
library(Matrix)


##### Part A) Deducing chemical formula #####

## known MW of elements
C_wt <- 12.0107
H_wt <- 1.00784
N_wt <- 14.0067
O_wt <- 15.999
P_wt <- 30.973762
S_wt <- 32.065
F_wt <- 18.998403
Cl_wt <- 35.453
Br_wt <- 79.904
Mol_name <- c("C","N","O","F","P","S","Cl","Br")
Mol_wt <- c(C_wt,N_wt,O_wt,F_wt,P_wt,S_wt,Cl_wt,Br_wt)
Mol_nElements <- length(Mol_wt)

## threshold for weight matching (e.g. 0.0001 means within 0.01% of target MW)
weight_threhold <- 0.0001


formula_deduction <- function(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount){

  ## setup of general variables
  Mol_count <- rep(0,Mol_nElements)
  Mol_count[0] <- HeavyAtomCount

  ## check naively a pure carbon compound does not exceed target MW
  if(sum(Mol_count*Mol_wt)>HeavyAtomMolWt){
    print("Houston, we have a problem. Even pure carbon is heavier than the given Molecular Weight")
    return(-1)
  } else {
    print("Safe for takeoff!")
  }

  ## matrix calculations to generate all possible chemical formulas
  composition <- t(compositions(HeavyAtomCount,Mol_nElements))
  weight_matrix <- composition %*% diag(Mol_wt)

  ## filter for formulas that match closely to target MW
  likely_composition <- which(rowSums(weight_matrix) <= (1+weight_threhold)*HeavyAtomMolWt & rowSums(weight_matrix) >= (1-weight_threhold)*HeavyAtomMolWt)
  candidates <- composition[likely_composition,]

  ## typically, we have more than one candidate within expected MW range
  ## the reason is that MW of N and O are very close to integers
  ## to narrow down the list, we make use of NO_count
  top_candidate <- which(candidates[,2]+candidates[,3]==NO_count)
  candidates[top_candidate,]
  if(length(top_candidate)>1){
    candidates <- candidates[top_candidate,]
    top_candidate <- which(candidates[,1]==HeavyAtomCount-HeteroatomCount)
  }
  rowSums(weight_matrix[likely_composition,])

  ## return result based on whether we found 0, 1, or >1 candidates
  if(length(top_candidate)==1){
    print("Congrats! There is only one plausible candidate composition:")
    correct_composition <- candidates[top_candidate,]
    names(correct_composition) <- c(Mol_name)
    return(correct_composition)

  } else if(length(top_candidate)>=1){
    print("Hmm... There are more than one plausible candidate composition, please ask a human for help:")
    candidate_list <- candidates[top_candidate,]
    colnames(candidate_list) <- Mol_name
    return(candidate_list)

  } else {
    print("Uh oh! No plausible candidates. Consider changing the weight matching threshold or rare elements.")
    return(-1)
  }
}


# toy example with 3 heavy atoms, Ethanol: 2*C,1*O
HeavyAtomCount <- 3
HeavyAtomMolWt <- 40.0204
NO_count <- 1
HeteroatomCount <- 1
formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)

# example from compound set 1, id=1
HeavyAtomCount <- 30
HeavyAtomMolWt <-423.772
NO_count <- 6
HeteroatomCount <- 8
formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)

# example from compound set 1, id=2
HeavyAtomCount <- 23
HeavyAtomMolWt <- 312.269
NO_count <- 6
HeteroatomCount <- 7
formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)

# example from compound set 1, id=3
HeavyAtomCount <- 26
HeavyAtomMolWt <- 340.22
NO_count <- 4
HeteroatomCount <- 6
formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)

# example from compound set 1, id=20
HeavyAtomCount <- 36
HeavyAtomMolWt <- 503.325
NO_count <- 8
HeteroatomCount <- 12
formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)

# example from compound set 1, id=22
HeavyAtomCount <- 30
HeavyAtomMolWt <- 495.168
NO_count <- 8
HeteroatomCount <- 13
formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)

# example from compound set 1, id=118
HeavyAtomCount <- 37
HeavyAtomMolWt <- 464.359
NO_count <- 8
HeteroatomCount <- 8
formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)

compound_set_1 <- read_excel("~/Downloads/Hacks/Challenge 1/Challeng_1_data/data/compound_set_1.xlsx")

compound_deducer <- function(id){
  row <- id+1
  data <- compound_set_1[row,]
  HeavyAtomCount <- data[[5]]
  HeavyAtomMolWt <- data[[6]]
  NO_count <- data[[10]]
  HeteroatomCount <- data[[19]]

  ## 40 atoms takes a while to compute
  ## skip larger molecules or risk memory overflow
  if(HeavyAtomCount>39){
    return(-1)
  }

  res <- formula_deduction(HeavyAtomMolWt,HeavyAtomCount,NO_count,HeteroatomCount)
  print(length(res))
  if(length(res)==Mol_nElements){
    return(1)
  } else if(length(res)>Mol_nElements){
    return(0)
  } else {
    return(-1)
  }
}


successes <- 0
almosts <- 0
fails <- 0
for(id in 1:nrow(compound_set_1)){
  outcome <- compound_deducer(id)
  if(outcome==1){
    successes <- successes+1
  } else if(outcome==0){
    almosts <- almosts+1
  } else {
    fails <- fails+1
  }
}

successes
almosts
fails


##### Part A) Summary #####

## we can successfull

##### Part B) Constructing adjacency matrices #####





# just toy example with carboxylic benzoic acid
correct_composition <- c(7,0,2,0,0,0,0)
H_available <- 6

# another toy example with toluene
#correct_composition <- c(7,0,0,0,0,0,0)
#H_available <- 8

# another toy example with di-carboxylic benzoic acid
#correct_composition <- c(8,0,4,0,0,0,0)
#H_available <- 6

# another toy example with phenol
#correct_composition <- c(6,0,1,0,0,0,0)
#H_available <- 6

element_list <- c()
for(element in 1:length(correct_composition)){
  print(element)
  element_list <- c(element_list,rep(Mol_name[element],correct_composition[element]))

}
element_list


benzene <- c()
benzene$mat <- read_excel("~/Downloads/benzene_mat.xlsx", na = "NA")
benzene$nAtoms <- 6
benzene$count <- 1
benzene$twoD <- c("C","C","C","C","C","C")

COOH <- c()
COOH$mat <- read_excel("~/Downloads/COOH_mat.xlsx", na = "NA")
COOH$nAtoms <- 3
COOH$count <- 1
COOH$twoD <- c("C","O","O")

phenol <- c()
phenol$mat <- read_excel("~/Downloads/phenol_mat.xlsx", na = "NA")
phenol$nAtoms <- 7
phenol$count <- 1
phenol$twoD <- c("C","C","C","C","C","C","O")



feature_names <- c("benzene","COOH","phenol")
feature_count <- c(benzene$count,COOH$count,phenol$count)
feature_nAtoms <- c(benzene$nAtoms,COOH$nAtoms,phenol$nAtoms)
feature_mat <- list(benzene$mat,COOH$mat,phenol$mat)
feature_twoD <- list(benzene$twoD,COOH$twoD,phenol$twoD)

feature_mapping <- c(
  rep(rep(feature_names[1],feature_nAtoms[1]),feature_count[1]),
  rep(rep(feature_names[2],feature_nAtoms[2]),feature_count[2]),
  rep(rep(feature_names[3],feature_nAtoms[3]),feature_count[3])
)

element_list <- c()
big_matrix <- c()
diag_idx <- 0
for(moiety in 1:length(feature_count)){
  left_padding <- diag_idx
  right_padding <- sum(correct_composition)-(diag_idx+feature_nAtoms[[moiety]])
  print(left_padding)
  print(right_padding)
  #print(moiety)
  new_block <- cbind(matrix(NA,feature_nAtoms[[moiety]],left_padding),as.matrix(feature_mat[[moiety]]),matrix(NA,feature_nAtoms[[moiety]],right_padding))
  big_matrix <- rbind(big_matrix,new_block)
  diag_idx <- diag_idx+feature_nAtoms[[moiety]]
  element_list <- c(element_list,feature_twoD[[moiety]])
}

big_matrix
colnames(big_matrix) <- feature_mapping
rownames(big_matrix) <- feature_mapping
element_list

for(row in 1:nrow(big_matrix)){
  if(element_list[row]=="O"){
    if(sum(big_matrix[row,],na.rm=TRUE)==2){
      big_matrix[row,][which(is.na(big_matrix[row,]))] <- 0
      big_matrix[,row][which(is.na(big_matrix[,row]))] <- 0
    } else if(rownames(big_matrix)[row]=="COOH"||rownames(big_matrix)[row]=="phenol"){
      if(sum(big_matrix[row,],na.rm=TRUE)==1){
        big_matrix[row,][which(is.na(big_matrix[row,]))] <- 0
        big_matrix[,row][which(is.na(big_matrix[,row]))] <- 0
        print(rownames(big_matrix)[row])
        H_available <- H_available-1
      }
    }
  }
  if(element_list[row]=="C"){
    if(sum(big_matrix[row,],na.rm=TRUE)==4){
      big_matrix[row,][which(is.na(big_matrix[row,]))] <- 0
      big_matrix[,row][which(is.na(big_matrix[,row]))] <- 0
    }
  }
  if(element_list[row]=="N"){
    if(sum(big_matrix[row,],na.rm=TRUE)==3){
      big_matrix[row,][which(is.na(big_matrix[row,]))] <- 0
      big_matrix[,row][which(is.na(big_matrix[,row]))] <- 0

    }
  }
}



big_matrix
H_available

# attempt to fill in missing spots by backtracking

next_blank <- which(is.na(big_matrix))[1]
next_col<- next_blank%%9
next_row <- floor(next_blank/9)

lower <- tril(big_matrix)
next_blank <- which(is.na(lower))

print(next_blank)
x_temp <- length(next_blank)
m_temp <- H_available

if(x_temp < m_temp){
  print("Uh oh! We don't have any place to put the Extra Hydrogens!")
} else {
  print("Testing the best placements for the Extra Hydrogens now...")
}

test_placements <- t(combn(x_temp,m_temp,function(x) replace(numeric(x_temp),x,1)))
placements <- 1 - test_placements

possible_matrices <- list()
for(row in 1:nrow(placements)){
  test_matrix <- big_matrix
  test_matrix[next_blank]<-placements[row,]

  full_matrix <- Matrix::forceSymmetric(test_matrix,uplo="L")
  #print(full_matrix)
  possible_matrices <- append(possible_matrices,list(full_matrix))

}

atomic_number_list <- c(6,6,6,6,6,6,6,8,8)


as.matrix(possible_matrices[[1]])

write.table(as.matrix(possible_matrices[[1]]),"~/Downloads/possible_matrix1.csv",sep=",",row.names = F, col.names=F)
write.table(atomic_number_list,"~/Downloads/possible_matrix1_headers.csv",sep=",",row.names = F, col.names=F)





# function for molecular feature contrstruction into adjacency matrix
benzene <- c()
benzene$mat <- read_excel("~/Downloads/benzene_mat.xlsx", na = "NA")
benzene$nAtoms <- 6


# function to use


#####











match_list <- list()

for(row_compos in c(1:nrow(composition))){
  Mol_count <- composition[row_compos,]
  #print(Mol_count)
  Total_wt <- sum(Mol_count*Mol_wt)
  if(abs(1-(Total_wt/HeavyAtomMolWt))<0.01){
    print("Hallelujah! Match is found")
    print(Mol_count*Mol_wt)
    print(Total_wt)
    match_list <- c(match_list,list(Mol_count))
  }
}
match_list










#combos(length(Mol_count))$binary[]

#temp <- combos(2)$ragged

#combo_matrix <- rbind()

matchfound <- c()
for(n_carbons in c(C_count:1)){
  for(row_combos in c(1:nrow(combos(HeavyAtomCount-1)$binary[]))){
    Mol_count <- c(n_carbons,combos(HeavyAtomCount-1)$binary[row_combos,])
    #print(Mol_count)
    Total_wt <- sum(Mol_count*Mol_wt)
    if(abs(1-(Total_wt/HeavyAtomMolWt))<0.01){
      matchfound <- Mol_count
      print("Hallelujah! Match is found")
    }
  }

}

print(matchfound)



while(!matchfound){
  Total_wt <- sum(Mol_count*Mol_wt)
  if(abs(1-(Total_wt/HeavyAtomMolWt))<0.01){
    matchfound <- TRUE
  }
  else {

  }
}
