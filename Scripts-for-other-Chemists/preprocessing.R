# Author: Haozhen Wu. 
# Date: 12/31/2015. Modified: 1/11/2015
# Convert data into required format.
# command line : Rscript preprocessing.R XX.csv(csv file name) 500(volume amount) 1(which version)

# Check if the system has dplyr package. If the system does not have the package, it will download by itself.
if("dplyr" %in% rownames(installed.packages())==FALSE){
  install.packages("dplyr")
  library(dplyr)
} else{
  library(dplyr)
}

# Based on each unique id of Synonyms, find the corresponding instances(row).
convert = function(i,volume, version){

# Find all the instances(row) that has the required volume amount of the Synonyms.
temp_500 = filter(alldata,Synonyms==i) %>% 
           filter( CA.Dose.Response.Data..Volume..nL.== volume |EC.Dose.Response.Data..Volume..nL.== volume|
                   PA.Dose.Response.Data..Volume..nL.== volume |SA.Dose.Response.Data..Volume..nL.== volume)

# Convert CA EC PA and SA of same Synonyms and required volume into same row.
result = cbind(
select(temp_500,Molecule.Name:CA.Dose.Response.Data....negative.control....) %>% 
  filter(CA.Dose.Response.Data..Volume..nL.== volume),
select(temp_500,EC.Dose.Response.Data..Run.Date:EC.Dose.Response.Data....negative.control....) %>%
  filter(EC.Dose.Response.Data..Volume..nL.== volume),
select(temp_500,PA.Dose.Response.Data..Run.Date:PA.Dose.Response.Data....negative.control....) %>%
  filter(PA.Dose.Response.Data..Volume..nL.== volume),
select(temp_500,SA.Dose.Response.Data..Run.Date:SA.Dose.Response.Data....negative.control....) %>%
  filter(SA.Dose.Response.Data..Volume..nL.== volume)
)

# Some Synonyms have different versions. Generate the output based on user specified version.
if (version == 1){
 result = result[1,] 
} else if (version ==2){
 if (dim(result)[1] == 2 ){
  result = result[2,] 
 }
}

return(result)
}

# Take command line argument.
args <- commandArgs(trailingOnly = TRUE)
alldata <- read.csv(args[1], stringsAsFactors=FALSE)
volumn = args[2]
version = args[3]
file_name <- unlist(strsplit(args[1], "[.]"))

# Find all the unique name of Synonyms.
syn = unique(alldata$Synonyms)
# Ignore the instances that do not have Synonyms names.
syn = syn[which(syn!="")]

# based on user specifed version, generate the corresponding result.
if (version == 1){
ver1 = sapply(syn, function(x) convert(x,volumn,version))
ver1 = t(ver1)
ver1 = data.frame(ver1)
ver1 <- data.frame(lapply(ver1, as.character), stringsAsFactors=FALSE)
#write.csv(ver1, file = paste0(file_name[1], "_version1.csv"))
write.csv(ver1, file = paste0(file_name[1],"_V",volumn, "_version1.csv"))

} else if(version ==2){
  ver2 = sapply(syn, function(x) convert(x,volumn,version))
  ver2 = t(ver2)
  ver2 = data.frame(ver2)
  ver2 <- data.frame(lapply(ver2, as.character), stringsAsFactors=FALSE)
  #write.csv(ver2, file = paste0(file_name[1], "_version2.csv"))
  write.csv(ver2, file = paste0(file_name[1],"_V",volumn, "_version2.csv"))
}



