#Author: Phoenix Logan, Haozhen Wu. 
# Usage: Rscript .csv 500(vol)


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
    return(y)
  }
  else{
    return(0)
  }
} 

#Plot the 4 dataset in Landscape view
plot_landscape <- function(df_name, Name){
  plot(df_name$neg_control, ylab= "", xlab= "", xaxt='n', 
       main= Name, type = "l", col = "black", ylim=c(0,100), cex=1.0, lwd=2)
  mtext("% Inhibition", side=2, line=2,las=0, font=1)
  abline(h=50,col="blue")# revised by haozhen,12/11
  par(new=TRUE)
  plot(df_name$mass, type="l", col= "red", xlab= NA, ylab=NA, axes=F, 
       lwd=1, ylim=c(0,0.5))
  axis(side=4, col.ticks="red", col.axis="red", las=0)
  axis(side=1, at=c(1:80), labels=ind$Wells, las=2, cex.axis=0.5)
  mtext("Weight by ESLD (mg)", las=0, side=4,line=2, col="red", cex=0.8) 
}

#Plot all four strains overlayed
plot_overlay <- function(CA,EC,PA,SA, unique_strain){
  #function will overlay all lines crated through % inhibition
  png(paste0(unique_strain,"_V",vol, "_overlay.png"), units = "in", width=11, height=8.5, res=300)
  par(mfrow=c(1,1), mar=c(2.0, 4.0,2.0,4.0), oma=c(0,0,3,1))
  title <- paste0(unique_strain, " at ",vol," nL")
  plot(CA$neg_control, ylab= "% Inhibition", xlab= "", xaxt='n', main=title, 
       type = "l", col = "purple", ylim=c(0,100), cex=1.0, lwd=4)
  abline(h=50,col="black") # revised by haozhen,12/11
  par(new=TRUE)
  lines(EC$neg_control, ylab= "% Inhibition", xlab= "", xaxt='n',main= "EC Sample",                  
        type = "l", col = "green", ylim=c(0,100), lwd=4)
  par(new=TRUE)
  lines(PA$neg_control, ylab= "", xlab= "",xaxt='n', main= "PA Sample",                  
        type = "l", col = "blue", ylim=c(0,100), lwd=4)
  par(new=TRUE)
  lines(SA$neg_control, ylab= "", xlab= "", xaxt='n',main= "SA Sample",                  
        type = "l", col = "orange", ylim=c(0,100), lwd=4)
  par(new=TRUE)
  plot(PA$mass, type="l", col= "red", xlab= NA, ylab=NA, axes=F, lwd=4, ylim=c(0,0.8))
  axis(side=4, col.ticks="red", col.axis="red", las=0)
  axis(side=1, at=c(1:80), labels=ind$Wells, font=0.2, las=2, cex.axis=0.5)
  mtext("Weight by ESLD", las=0, side=4,line=3, col="red", cex=0.8)
  # revised by Haozhen Wu. 12/4.
  legend("topright", c("CA","EC","PA","SA", "Mass","50%"), col=c("purple", "green", "blue","orange", "red","black"), cex=1, pch=19, pt.cex =0.5)
  dev.off()
}

#Parse data by each unique strain and plot landscape and overlay graphs
get_strain_df <- function(unique_strain, prefix, wells, strain){
  strain_ID <- paste0(prefix, ".Dose.Response.Data....negative.control....")
  strain_volume <- paste0(prefix, ".Dose.Response.Data..Volume..nL.")
  neg <- unique_strain[,strain_ID]
  inv_neg <- as.vector(sapply(neg, inv_neg))
  #new_df <- data.frame(Wells=wells, mass=as.numeric(unique_strain$Mass..mg.), volume=as.numeric(unique_strain[,strain_volume]), neg_control=as.numeric(inv_neg), stringsAsFactors = FALSE)
  # revised by haozhen 16/1/13
  new_df <- data.frame(Molecule_Name=unique_strain$Molecule.Name,Wells=wells, mass=as.numeric(unique_strain$Mass..mg.), volume=as.numeric(unique_strain[,strain_volume]), neg_control=as.numeric(inv_neg), stringsAsFactors = FALSE)
  
  data_500 <- which(new_df$volume == vol)
  DATA_parsed <- new_df[data_500,]
  # revised by haozhen 12/30
  if(length(which(DATA_parsed$Wells %in% ind$Wells)) !=80){
    complete_df <- merge(ind_2, DATA_parsed, by="Wells", sort=FALSE)
  } else {
    complete_df <- merge(ind, DATA_parsed, by="Wells", sort=FALSE)
  }
  #
  complete_df_larger50 = complete_df[which(complete_df$neg_control>=50),]
  
  # haozhen
  complete_df_larger50$Strain = rep(strain,dim(complete_df_larger50)[1])
  complete_df_larger50$index = rownames(complete_df_larger50)
  complete_df_larger50 = complete_df_larger50[,c(2,6,1,7,3,4,5)]
  #
  if(dim( complete_df_larger50)[1]>0){
    write.csv(complete_df_larger50, file = paste0(strain,"_",prefix,"_above50","_V",vol,".csv"),row.names = FALSE)  # revised by haozhen,12/11
    # revised by haozhen 16/01/12
    if(prefix == "CA"){
      write.table(complete_df_larger50, file = paste0("CA_above50_V",vol, ".csv"),append = T, row.names = FALSE, col.names = FALSE, sep = ",")
    } else if( prefix == "EC"){
      write.table(complete_df_larger50, file = paste0("EC_above50_V",vol, ".csv"),append = T, row.names = FALSE, col.names = FALSE,sep = ",")
    } else if (prefix == "PA"){
      write.table(complete_df_larger50, file = paste0("PA_above50_V",vol, ".csv"),append = T, row.names = FALSE, col.names = FALSE,sep = ",")
    } else if( prefix == "SA"){
      write.table(complete_df_larger50, file = paste0("SA_above50_V",vol, ".csv"),append = T, row.names = FALSE, col.names = FALSE,sep = ",")
    }
  }
  #
  write.csv(complete_df, file=paste0(strain,"_", prefix,"_V",vol, ".csv"))
  return(complete_df)
}

get_data <- function(unique_strain, df){
  #grep out the unique strain
  #unique_selection_1 <- grep(paste0(unique_strain, "-[A-Z][0-9][0-9]"), df$Synonyms)
  unique_selection_1 <- grep(paste0(unique_strain, "-[A-Z]*"), df$Synonyms) # revised by haozhen 12/30
  unique_2 <- df[unique_selection_1,]
  unique_syns <- unique_2$Synonyms
  wells <- as.vector(sapply(unique_syns, get_strain_info, num=2))
  #CA Data
  CA_2 <- get_strain_df(unique_2, "CA", wells, unique_strain) 
  #EC Data
  EC_2 <- get_strain_df(unique_2, "EC", wells, unique_strain)
  #PA Data
  PA_2 <- get_strain_df(unique_2, "PA", wells, unique_strain)
  #SA Data
  SA_2 <- get_strain_df(unique_2, "SA", wells, unique_strain)
  
  png(file=paste0(unique_strain,"_V",vol, ".png"), units = "in", width=11, height=8.5, res=300)
  par(mfrow=c(4,1), mar=c(2.0, 4.0,2.0,4.0), oma=c(1,1,3,1))
  plot_landscape(CA_2, "CA")
  plot_landscape(EC_2, "EC")
  plot_landscape(PA_2, "PA")
  plot_landscape(SA_2, "SA")
  mtext(paste0(unique_strain,"_V",vol," nL"), side=3, line=1, outer=TRUE, cex=2, font=2)
  dev.off()
  plot_overlay(CA_2,EC_2,PA_2,SA_2, unique_strain)
}


#create interactive command line arguments later
args <- commandArgs(trailingOnly = TRUE)
data_table <- read.csv(args[1], stringsAsFactors=FALSE)
vol = args[2]
#to get well names 
synonyms <- data_table$Synonyms

#get file name
file_name <- unlist(strsplit(args[1], "[.]"))
jpeg(paste0(file_name[1], "landscape.jpg"))

#get strain name and well
well_idx_vector <- as.vector(sapply(synonyms, get_strain_info, num=2))
well_strain_vector <- as.vector(sapply(synonyms, get_strain_info, num=1))

#Organize the data by ascending/descending order
ind <- data.frame(Wells=paste0(rep(LETTERS[1:8], each=10), c("02","03","04","05","06","07","08","09","10","11", "11","10","09","08","07","06","05","04","03","02")), stringsAsFactors=FALSE)
ind_2<-data.frame(Wells=paste0(rep(LETTERS[1:8], each=10), c("2","3","4","5","6","7","8","9","10","11", "11","10","9","8","7","6","5","4","3","2")), stringsAsFactors=FALSE)
unique_strains <- unique(well_strain_vector)
# revised by haozhen 16/01/12
header = data.frame(Molecule_Name = character(0),Strain = character(0), Wells = character(0), index = character(0) , mass = numeric(0), volume = numeric(0), neg_control = numeric(0))
write.csv(header,file = paste0("CA_above50_V",vol, ".csv"),row.names = F)
write.csv(header,file = paste0("EC_above50_V",vol, ".csv"),row.names = F)
write.csv(header,file = paste0("PA_above50_V",vol, ".csv"),row.names = F)
write.csv(header,file = paste0("SA_above50_V",vol, ".csv"),row.names = F)
#
sapply(unique_strains, get_data, df=data_table)  


paste0("CA_above50_V",vol, ".csv")
