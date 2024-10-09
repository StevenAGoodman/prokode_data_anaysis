# filter empty stuff
# done in powershell
library(GEOquery)
library(stringr)
library(fixr)
library(TAF)

rmdir("C:/Users/steve/PROKODE-DOCKER/prokode/beta_training/GEO_expression_data", recursive = T)

# check each for negative numbers or wierd number ranges
file_names <- list.files(path = "C:/Users/steve/PROKODE-DOCKER/prokode/beta_training/GEO_expression_data", pattern = "*.csv", full.names = TRUE, recursive = T)

for (file in file_names){
    gse_id = str_extract(file, "(GSE\\d+)", group = 1)
    print(gse_id)
    gse = getGEO(gse_id, GSEMatrix =F)
    gsm = GSMList(gse)[[1]]
    organism = Meta(gsm)$organism_ch1
    organism = word(organism, 1,2, sep=" ")
    print(organism)

    dir.create(file.path(substr(file, 1, 71), organism))

    outp = paste0(substr(file, 1, 72), organism, '/', gse_id, ".csv")
    print(outp)
    file.rename(from = file, to = outp) 
}


#Loop through each file
for (file in file_names) {
    data_df = read.csv(file, header = T)
    data_df = data_df[2:length(data_df)]
    val_range = range(data_df, na.rm = T)
    print(val_range)
    if (length(check_for_negative_values(data_df)) != 0){
        file.rename(from = file, to = paste0(substr(file, 1, 85), "/bad", substr(file,86,nchar(file))))
    }
    if (val_range[2] < 25){
        file.rename(from = file, to = paste0(substr(file, 1, 85), "/bad", substr(file,86,nchar(file))))     
    }
}

file_names =list.files(path = "C:/Users/steve/PROKODE-DOCKER/prokode/beta_training/GEO_expression_data", pattern = "*.log", full.names = TRUE, recursive = T)

for (file in file_names){
    file.rename(from = file, to = paste0(substr(file, 1, 72), "/logs", substr(file,72,nchar(file))))     
}
    # if there are, verify that its a ration theyre measuring and attempt to find raw intensities

# organize columns into groups of same experimental condition
    # rename ones without "| n" to have "| 0"
# organize into organism 