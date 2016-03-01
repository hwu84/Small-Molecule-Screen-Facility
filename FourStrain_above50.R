# Author: Haozhen Wu
# Date: 1/28/2016. Modified: 2/2/2016
# input: xx(CA,EC,PA,SA)_above50_V(50,100,250,500).csv 
#find molecule overlap among 4 strains(CA,EC,PA,SA) and sort them in decreasing order


get_overlap = function(x){                                                                                          
  k = 1
  for (i in molecule_name){
    # check molecule satisfy all 4 strains
    if( i %in% CA & i %in% EC & i %in% PA & i %in% SA){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = SA_above50$Strain[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = SA_above50$Wells[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = SA_above50$index[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = SA_above50$mass[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "yes"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "yes"
      four_strains_above50$SA[k] = "yes"
      k = k + 1
    }
    
  }
  #molecule_name = molecule_name[!four_strains_above50$Molecule_Name %in% molecule_name]
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy CA EC PA
    if( i %in% CA & i %in% EC & i %in% PA){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = EC_above50$Strain[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = EC_above50$Wells[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = EC_above50$index[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = EC_above50$mass[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "yes"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "yes"
      four_strains_above50$SA[k] = "no"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy CA EC SA
    if( i %in% CA & i %in% EC & i %in% SA){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = SA_above50$Strain[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = SA_above50$Wells[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = SA_above50$index[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = SA_above50$mass[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "yes"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "no"
      four_strains_above50$SA[k] = "yes"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy  EC PA SA
    if( i %in% EC & i %in% PA & i %in% SA){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = SA_above50$Strain[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = SA_above50$Wells[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = SA_above50$index[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = SA_above50$mass[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "no"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "yes"
      four_strains_above50$SA[k] = "yes"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy CA EC
    if( i %in% CA & i %in% EC ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = CA_above50$Strain[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = CA_above50$Wells[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = CA_above50$index[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = CA_above50$mass[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "yes"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "no"
      four_strains_above50$SA[k] = "no"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy CA PA
    if( i %in% CA & i %in% PA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = PA_above50$Strain[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = PA_above50$Wells[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = PA_above50$index[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = PA_above50$mass[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "yes"
      four_strains_above50$EC[k] = "no"
      four_strains_above50$PA[k] = "yes"
      four_strains_above50$SA[k] = "no"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy CA SA
    if( i %in% CA & i %in% SA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = SA_above50$Strain[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = SA_above50$Wells[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = SA_above50$index[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = SA_above50$mass[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "yes"
      four_strains_above50$EC[k] = "no"
      four_strains_above50$PA[k] = "no"
      four_strains_above50$SA[k] = "yes"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy EC PA
    if( i %in% EC & i %in% PA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = PA_above50$Strain[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = PA_above50$Wells[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = PA_above50$index[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = PA_above50$mass[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "no"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "yes"
      four_strains_above50$SA[k] = "no"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  four_strains_above50$Molecule_Name[1:(k-1)]
  for (i in molecule_name){
    # check molecule satisfy EC SA
    if( i %in% EC & i %in% SA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = SA_above50$Strain[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = SA_above50$Wells[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = SA_above50$index[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = SA_above50$mass[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "no"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "no"
      four_strains_above50$SA[k] = "yes"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy PA SA
    if( i %in% PA & i %in% SA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = SA_above50$Strain[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = SA_above50$Wells[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = SA_above50$index[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = SA_above50$mass[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "no"
      four_strains_above50$EC[k] = "no"
      four_strains_above50$PA[k] = "yes"
      four_strains_above50$SA[k] = "yes"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy CA
    if( i %in% CA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = CA_above50$Strain[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = CA_above50$Wells[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = CA_above50$index[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = CA_above50$mass[which(CA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "yes"
      four_strains_above50$EC[k] = "no"
      four_strains_above50$PA[k] = "no"
      four_strains_above50$SA[k] = "no"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy EC
    if( i %in% EC ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = EC_above50$Strain[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = EC_above50$Wells[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = EC_above50$index[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = EC_above50$mass[which(EC_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "no"
      four_strains_above50$EC[k] = "yes"
      four_strains_above50$PA[k] = "no"
      four_strains_above50$SA[k] = "no"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy PA
    if( i %in% PA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = PA_above50$Strain[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = PA_above50$Wells[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = PA_above50$index[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = PA_above50$mass[which(PA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "no"
      four_strains_above50$EC[k] = "no"
      four_strains_above50$PA[k] = "yes"
      four_strains_above50$SA[k] = "no"
      k = k + 1
    }
    
  }
  molecule_name = molecule_name[!(molecule_name %in% four_strains_above50$Molecule_Name)]
  
  for (i in molecule_name){
    # check molecule satisfy SA
    if( i %in% SA ){
      four_strains_above50$Molecule_Name[k] = i
      four_strains_above50$Strain[k] = SA_above50$Strain[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$Wells[k] = SA_above50$Wells[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$index[k] = SA_above50$index[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$mass[k] = SA_above50$mass[which(SA_above50$Molecule_Name == i)]
      four_strains_above50$CA[k] = "no"
      four_strains_above50$EC[k] = "no"
      four_strains_above50$PA[k] = "no"
      four_strains_above50$SA[k] = "yes"
      k = k + 1
    }
    
  }
  return(four_strains_above50)
}


## input csv files should follow CA,EC,PA,SA
#CA_above50 <- read.csv("~/Desktop/study/SmallMolecular/R/Rscript/GENE/sort_u19_4strains_all_data_12_30_2015_GA/sort_ver1/V50/CA_above50_V50.csv",stringsAsFactors = F)
#EC_above50 <- read.csv("~/Desktop/study/SmallMolecular/R/Rscript/GENE/sort_u19_4strains_all_data_12_30_2015_GA/sort_ver1/V50/EC_above50_V50.csv",stringsAsFactors = F)
#PA_above50 <- read.csv("~/Desktop/study/SmallMolecular/R/Rscript/GENE/sort_u19_4strains_all_data_12_30_2015_GA/sort_ver1/V50/PA_above50_V50.csv",stringsAsFactors = F)
#SA_above50 <- read.csv("~/Desktop/study/SmallMolecular/R/Rscript/GENE/sort_u19_4strains_all_data_12_30_2015_GA/sort_ver1/V50/SA_above50_V50.csv",stringsAsFactors = F)

args <- commandArgs(trailingOnly = TRUE)
#data_table <- read.csv(args[1], stringsAsFactors=FALSE)

CA_above50 = read.csv(args[1], stringsAsFactors=FALSE)
EC_above50 = read.csv(args[2], stringsAsFactors=FALSE)
PA_above50 = read.csv(args[3], stringsAsFactors=FALSE)
SA_above50 = read.csv(args[4], stringsAsFactors=FALSE)
vol = args[5]

# Get molecule name of each strain
CA <- CA_above50$Molecule_Name
EC <- EC_above50$Molecule_Name
PA <- PA_above50$Molecule_Name
SA <- SA_above50$Molecule_Name

# get union of molecule from 4 strains
molecule_name = Reduce(union, list(CA,EC,PA,SA))
# create dataframe.
four_strains_above50 = data.frame(Molecule_Name = character(length(molecule_name)), 
                                  Strain = character(length(molecule_name)),
                                  Wells = character(length(molecule_name)),
                                  index = numeric(length(molecule_name)),
                                  mass = numeric(length(molecule_name)),
                                  CA = character(length(molecule_name)),
                                  EC = character(length(molecule_name)),PA = character(length(molecule_name)),
                                  SA = character(length(molecule_name)),stringsAsFactors = FALSE)

result = get_overlap(molecule_name)
write.csv(result,file = paste0("FourStrains_above50_V",vol, ".csv"),row.names = F)


