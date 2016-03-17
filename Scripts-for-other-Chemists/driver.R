#Author: Haozhen Wu. 3/11
# Usage: Rscript driver.R inputfile(eg. export_mar11.csv
#
#

#strip off the first part of the Synonym for well indexes
get_strain_info <- function(syn, num){
  well_idx_str <- strsplit(syn, "-")
  return(well_idx_str[[1]][num])
}

#get % Inhibition
inv_neg <- function(x) {
  y <- 100-x
  #remove elements less than zero
  if (y > 0){
    return(round(y,2))
  }
  else{
    return(0)
  }
} 

#Plot the individual plot, eg. PA at vol=500
plot_landscape <- function(df_name, Name,color){
  plot(df_name$neg_control, ylab= "", xlab= "", xaxt='n', 
       main= Name, type = "l", col = color, ylim=c(0,100), cex=1.0, lwd=2)
  mtext("% Inhibition", side=2, line=2,las=0, font=1,col = color)
  abline(h=50,col="black")
   axis(side=1, at=c(1:80), labels=ind$Wells, las=2, cex.axis=0.5)
  }
# plot the mass 
plot_mass <- function(df_name,Name,color){
  plot(df_name$ELSD, type="l", col= color, xlab= "", ylab="", xaxt='n',
            main = Name,lwd=1, ylim=c(0,0.5)) 
  axis(side=1, at=c(1:80), labels=ind$Wells, las=2, cex.axis=0.5) 
  mtext("Weight by ESLD (mg)", las=0, side=2,line=2, col="black", cex=0.8)  

}

#Plot all four strains overlayed
plot_overlay <- function(CA,EC,PA,SA, unique_strain,vol){
  #function will overlay all lines crated through % inhibition
  png(paste0("./vol",vol,"/",unique_strain,"_V",vol, "_overlay.png"), units = "in", width=11, height=8.5, res=300)
  par(mfrow=c(1,1), mar=c(2.0, 4.0,2.0,4.0), oma=c(0,0,3,1))
  title <- paste0(unique_strain, " at ",vol," nL")
  plot(CA$neg_control, ylab= "% Inhibition", xlab= "", xaxt='n', main=title, 
       type = "l", col = "green", ylim=c(0,100), cex=1.0, lwd=4)
  abline(h=50,col="black") # revised by haozhen,12/11
  par(new=TRUE)
  lines(EC$neg_control, ylab= "% Inhibition", xlab= "", xaxt='n',main= "EC Sample",                  
        type = "l", col = "red", ylim=c(0,100), lwd=4)
  par(new=TRUE)
  lines(PA$neg_control, ylab= "", xlab= "",xaxt='n', main= "PA Sample",                  
        type = "l", col = "blue", ylim=c(0,100), lwd=4)
  par(new=TRUE)
  lines(SA$neg_control, ylab= "", xlab= "", xaxt='n',main= "SA Sample",                  
        type = "l", col = "orange", ylim=c(0,100), lwd=4)
  par(new=TRUE)
  plot(PA$ELSD, type="l", col= "black", xlab= NA, ylab=NA, axes=F, lwd=4, ylim=c(0,0.8))
  axis(side=4, col.ticks="red", col.axis="red", las=0)
  axis(side=1, at=c(1:80), labels=ind$Wells, font=0.2, las=2, cex.axis=0.5)
  mtext("Weight by ESLD", las=0, side=4,line=3, col="red", cex=0.8)
  # revised by Haozhen Wu. 12/4.
  legend("topright", c("PA","EC","CA","SA", "Mass","50%"), col=c("blue", "red", "green","orange", "black","black"), cex=1, pch=19, pt.cex =0.5)
  dev.off()
}

# Prepare dataset for plot for each strain
get_strain_df <- function(unique_strain, prefix, wells, strain,vol){
  strain_ID = paste0(prefix,"_",vol)
  neg <- unique_strain[,strain_ID]
  inv_neg <- as.vector(sapply(neg, inv_neg))
  new_df = data.frame(SMSSF_ID = unique_strain$SMSSF_ID, Wells = wells, ELSD = as.numeric(unique_strain$ELSD),volume = as.numeric(rep(vol,dim(unique_strain)[1])), neg_control = as.numeric(inv_neg), stringsAsFactors = F)
  
  DATA_parsed <- new_df
  # revised by haozhen 12/30
  if(length(which(DATA_parsed$Wells %in% ind$Wells)) !=80){
    complete_df <- merge(ind_2, DATA_parsed, by="Wells", sort=FALSE)
  } else {
    complete_df <- merge(ind, DATA_parsed, by="Wells", sort=FALSE)
  }
  
  return(complete_df)
}

# Driver method. It will call --> get_starin_df --> plot_landscape --> plot_mass --> plot_overlay 
get_data <- function(unique_strain, df,vol){
  #grep out the unique strain
  #unique_selection_1 <- grep(paste0(unique_strain, "-[A-Z][0-9][0-9]"), df$Synonyms)
  unique_selection_1 <- grep(paste0(unique_strain, "-[A-Z]*"), df$X96_Plate_Well) # revised by haozhen 12/30
  unique_2 <- df[unique_selection_1,]
  unique_syns <- unique_2$X96_Plate_Well
  wells <- as.vector(sapply(unique_syns, get_strain_info, num=2))
  #CA Data
  CA_2 <- get_strain_df(unique_2, "CA", wells, unique_strain,vol) 
  #EC Data
  EC_2 <- get_strain_df(unique_2, "EC", wells, unique_strain,vol)
  #PA Data
  PA_2 <- get_strain_df(unique_2, "PA", wells, unique_strain,vol)
  #SA Data
  SA_2 <- get_strain_df(unique_2, "SA", wells, unique_strain,vol)
  
  png(file=paste0("./vol",vol,"/",unique_strain,"_V",vol, ".png"), units = "in", width=11, height=8.5, res=300)
  par(mfrow=c(5,1), mar=c(2.0, 4.0,2.0,4.0), oma=c(1,1,3,1))
  plot_landscape(PA_2, "PA","blue")
  plot_landscape(EC_2, "EC","red")
  plot_landscape(CA_2, "CA","green")
  plot_landscape(SA_2, "SA","orange")
  plot_mass(CA_2,"Mass","black")
  mtext(paste0(unique_strain,"_V",vol," nL"), side=3, line=1, outer=TRUE, cex=2, font=2)
  dev.off()
  plot_overlay(CA_2,EC_2,PA_2,SA_2, unique_strain,vol)
}

# generate sorted file and row larger than 50. 
generateFile <- function(fileName){
  
  names = c("PA_500","SA_500","EC_500","CA_500","PA_250","SA_250","EC_250","CA_250",
            "PA_100","SA_100","EC_100","CA_100","PA_50","SA_50","EC_50","CA_50")
  # convert value into negative control
  for(i in names){
    data_table[,which(colnames(data_table) == i)] = sapply( data_table[,which(colnames(data_table) == i)], inv_neg)
  }
  #
  well_idx_vector <- as.vector(sapply(synonyms, get_strain_info, num=2))
  well_strain_vector <- as.vector(sapply(synonyms, get_strain_info, num=1))
  unique_strain = unique(well_strain_vector)
  data_table$strain = well_strain_vector
  data_table$Wells = well_idx_vector
  #file_name <- unlist(strsplit(args[1], "[.]"))
  
  boolen = F
  boolen2 = F
  for( i in unique_strain){
    temp_df = data_table[which(data_table$strain %in% i),] 
    if(length(which(temp_df$Wells %in% ind$Wells)) !=80){
      temp_df <- merge(ind_2, temp_df, by="Wells", sort=FALSE)
    } else {
      temp_df <- merge(ind, temp_df, by="Wells", sort=FALSE)
    }
    temp_df$strain = NULL
    temp_df$Wells = NULL
    
    index_larger50 = which(temp_df$PA_500 >=50 | temp_df$SA_500 >=50 |temp_df$EC_500 >=50 |temp_df$CA_500 >=50 |
                             temp_df$PA_250 >=50 | temp_df$SA_250 >=50 |temp_df$EC_250 >=50 |temp_df$CA_250 >=50 |
                             temp_df$PA_100 >=50 | temp_df$SA_100 >=50 |temp_df$EC_100 >=50 |temp_df$CA_100 >=50 |
                             temp_df$PA_50 >=50 | temp_df$SA_50 >=50 |temp_df$EC_50 >=50 |temp_df$CA_50 >=50 
    )
    if(length(index_larger50 >0)){
      temp_df_larger50 = temp_df[index_larger50,]
      if(boolen2 == F){
        write.table(temp_df_larger50, file = paste0("./File/",fileName,"_larger50.csv"),sep = ",",row.names = F,col.names = T)
        #write.table(temp_df_larger50, file = paste0("_larger50.csv"),sep = ",",row.names = F,col.names = T)
        
        boolen2 = T
      } else{
        write.table(temp_df_larger50,file = paste0("./File/",fileName,"_larger50.csv"),append = T,sep = ",",row.names = F,col.names = F)
        #write.table(temp_df_larger50,file = paste0("_larger50.csv"),append = T,sep = ",",row.names = F,col.names = F)
        
        } 
      
    }
    
    if(boolen == F){
      write.table(temp_df, file = paste0("./File/",fileName,"_sorted.csv"),sep = ",", row.names = F, col.names = T)
      boolen = T
    }else{
      write.table(temp_df, file = paste0("./File/",fileName,"_sorted.csv"), append = T, sep = ",",row.names = F,col.names = F)
      }
  
    }
  
}

# function to remove duplicate and keep the latest one.
# assume when duplicate occur, latest one always on top of the others
remove_dup = function(data_table){
unique_id = unique(data_table$SMSSF_ID)
removed_num = 0
for(i in unique_id){
  if(length(which(data_table$SMSSF_ID %in% i)) >1){
   removed_pos = which(data_table$SMSSF_ID %in% i)[-1]
   data_table = data_table[-removed_pos,]
   removed_num = length(removed_pos) + removed_num
  }
}
cat("Removed duplicate: ",removed_num,"\n")
return(data_table)
}

########################################### Above are self defined methods #############


#create interactive command line arguments later
args <- commandArgs(trailingOnly = TRUE)
data_table <- read.csv(args[1], stringsAsFactors=FALSE)
# remove duplicate
data_table = remove_dup(data_table)

# Previous Synonyms now is X96_Plate_Well, same meaning
synonyms <- data_table$X96_Plate_Well


#get file name
file_name <- unlist(strsplit(args[1], "[.]"))[1]


#get strain name and well
well_idx_vector <- as.vector(sapply(synonyms, get_strain_info, num=2))
well_strain_vector <- as.vector(sapply(synonyms, get_strain_info, num=1))

#Organize the data by ascending/descending order
ind <- data.frame(Wells=paste0(rep(LETTERS[1:8], each=10), c("02","03","04","05","06","07","08","09","10","11", "11","10","09","08","07","06","05","04","03","02")), stringsAsFactors=FALSE)
ind_2<-data.frame(Wells=paste0(rep(LETTERS[1:8], each=10), c("2","3","4","5","6","7","8","9","10","11", "11","10","9","8","7","6","5","4","3","2")), stringsAsFactors=FALSE)
unique_strains <- unique(well_strain_vector)

# call method to generate sorted file and row with PA,EC,CA,SA value larger than 50.
generateFile(fileName = file_name)

# For each volumns, generate its individual plots and overlay plot
for(VOL in c(500,250,100,50)){
sapply(unique_strains, get_data, df=data_table, vol = VOL)  
}

