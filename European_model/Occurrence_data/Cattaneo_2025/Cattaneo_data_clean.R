### LITERATURE REVIEW DATA CLEANING ###

# The disease data in this script was retrieved from a literature review on Dengue, Chikungunya and Zika local outbreaks in Europe. 
# conducted by Cattaneo et al., feb 2025,the Lancet regional health: https://www.thelancet.com/journals/lanepe/article/PIIS2666-7762(25)00023-7/fulltext

# Data extraction method: They performed a systematic review of the literature published from January 1, 2007, to January 31, 2024, reporting autochthonous cases 
# of dengue, chikungunya, and Zika detected in Europe. They searched MEDLINE, EMBASE, and ECDC reports.

# The supplementary info they provided lists all articles and outbreaks included in their literature review in a PDF format
# The pdf was converted to a csv file, the data was then cleaned and processed in R, associating outbreak locations to their corresponding NUTS3 regions (following script)

library(raster)
library(sp)
library(readxl)
library(tidyr)

# Load European NUTS3 shapefile (MOOD project Willy Wint)
nuts3 = crop(shapefile("../moodpolygonsjanmasks21/VectornetDATAforMOODjan21.shp"), extent(-10,33,34.5,72)) # NUT3 polygons shapefile
nuts3 = subset(nuts3, !CountryISO%in%c("MA","DZ","TN","MT","TR","CI","MD","UA","RU","FO","IS","GE","BY","IS","GL","FO","CY","SJ"))
nuts3_data = as.data.frame(matrix(nrow = dim(nuts3@data)[1], ncol = 3))
colnames(nuts3_data) = c("country", "location_name", "NUTS3")
nuts3_data$NUTS3 = nuts3@data[,"LocationCo"]
nuts3_data$country = nuts3@data[,"CountryNam"]
nuts3_data$location_name = nuts3@data$LocationNa


# Read xlxs file from literature review (previously converted from pdf)
path = "/Users/kylaserres/Dropbox/PhD/Research projects/ENM-arboviral-projections/Scripts_and_Data/Data/Presence_data/Literature review/" 
excel_sheets(paste0(path,"included_occurrences_Cattaneo.xlsx"))
# sheet1
data = read_excel(paste0(path,"included_occurrences_Cattaneo.xlsx"), sheet = "disease_data_1")
data = data[,c("Location of the study (country)", "Outbreak dealt with (disease and year)", "Location/area involved in the outbreak",
               "Total number of autochthonous cases confirmed","Authors of the study")]
colnames(data) = c("country", "disease","location", "case_number", "authors")
#split disease column into disease and outbreak year
data = separate(data, col = disease, into = c("disease", "year"), sep = " ")
data$disease = gsub(",", "", data$disease)

#sheet2
data2 = read_excel(paste0(path,"included_occurrences_Cattaneo.xlsx"), sheet = "disease_data_2")
data2 = data2[,c("Country", "Outbreak", "Location/area involved in the outbreak", "Total number of autochthonous cases confirmed")]
colnames(data2) = c("country", "disease","location", "case_number")
data2$authors = NA
#split disease column into disease and outbreak year
data2 = separate(data2, col = disease, into = c("disease", "year"), sep = " ")
data2$disease = gsub(",", "", data2$disease)

#combine both sheets
all_data = rbind(data, data2)

# curation
  # remove all characters that are not letters
unique_chars = unique(unlist(strsplit(all_data$location, "")))
print(unique_chars)
all_data$location=gsub("[\n()/,]", "", all_data$location) 
all_data = all_data[!is.na(all_data$disease), ] # remove rows for which no disease info is available
all_data = all_data[!all_data$location=="NA", ] # remove rows for which no location info is available
all_data = all_data[!all_data$location=="Funchal and neighbouring", ] # remove rows for Maderia (PT) (not in continental europe)
all_data = all_data[!all_data$location=="Catalunya", ] # remove rows for Catalunya (NUTS2 region)
all_data$location = gsub("Department", "", all_data$location) # remove department word 
all_data$location = gsub("department", "", all_data$location) # remove department word 
all_data$location = gsub("province", "", all_data$location) # remove province word 
all_data$location = gsub("Province", "", all_data$location) # remove province word 

# Correct department (NUTS3) names 
all_data$location = ifelse(grepl("Var", all_data$location), "Var", all_data$location)
all_data$location = ifelse(grepl("Gard", all_data$location), "Gard", all_data$location)
all_data$location = ifelse(grepl("Hérault", all_data$location), "Herault", all_data$location)
all_data$location = ifelse(grepl("Herault", all_data$location), "Herault", all_data$location)
all_data$location = ifelse(grepl("Paris", all_data$location), "Paris", all_data$location)
all_data$location = ifelse(grepl("Limeil-Brévannes", all_data$location), "Paris", all_data$location)
all_data$location = ifelse(grepl("Rome", all_data$location), "Roma", all_data$location)
all_data$location = ifelse(grepl("Roma", all_data$location), "Roma", all_data$location)
all_data$location = ifelse(grepl("Anzio", all_data$location), "Roma", all_data$location)
all_data$location = ifelse(grepl("Murcia", all_data$location), "Murcia", all_data$location)
all_data$location = ifelse(grepl("Drôme", all_data$location), "Drome", all_data$location)
all_data$location = ifelse(grepl("Corsica", all_data$location), "Corse-du-Sud", all_data$location)
all_data$location = ifelse(grepl("Ravenna", all_data$location), "Ravenna", all_data$location)
all_data$location = ifelse(grepl("Bouches- du- Rhône", all_data$location), "Bouches-du-Rhone", all_data$location)
all_data$location = ifelse(grepl("Bouches-du-Rhône", all_data$location), "Bouches-du-Rhone", all_data$location)
all_data$location = ifelse(grepl("Bouches- du-Rhône", all_data$location), "Bouches-du-Rhone", all_data$location)
all_data$location = ifelse(grepl("Rhône-Alpes", all_data$location), "Rhone", all_data$location)
all_data$location = ifelse(grepl("Alpes", all_data$location), "Alpes-Maritimes", all_data$location)
all_data$location = ifelse(grepl("Lodi", all_data$location), "Lodi", all_data$location)
all_data$location = ifelse(grepl("Pyrénées", all_data$location), "Pyrenees-Orientales", all_data$location)
all_data$location = ifelse(grepl("Bigorre", all_data$location), "Hautes-Pyrenees", all_data$location)
all_data$location = ifelse(grepl("Gaude", all_data$location), "Alpes-Maritimes", all_data$location)
all_data$location = ifelse(grepl("Gattières", all_data$location), "Alpes-Maritimes", all_data$location)
all_data$location = ifelse(grepl("Montpellier", all_data$location), "Herault", all_data$location)
all_data$location = ifelse(grepl("Nîmes|Nimes", all_data$location), "Gard", all_data$location)
all_data$location = ifelse(grepl("Aix-en-Provence area", all_data$location), "Bouches-du-Rhone", all_data$location)
all_data$location = ifelse(grepl("Nice|Hyères|Fréjus|Frejus", all_data$location), "Alpes-Maritimes", all_data$location)
all_data$location = ifelse(grepl("La Croix-Valmer|Fayence", all_data$location), "Var", all_data$location)
all_data$location = ifelse(grepl("Perpignan", all_data$location), "Pyrenees-Orientales", all_data$location)
all_data$location = ifelse(grepl("Barcelona|Badalona", all_data$location), "Barcelona", all_data$location)
all_data$location = ifelse(grepl("Ibiza", all_data$location), "Eivissa y Formentera", all_data$location)
all_data$location = ifelse(grepl("Montecchio Maggiore", all_data$location), "Treviso", all_data$location)
all_data$location = ifelse(grepl("Tolosa|Montauban|La Salvetat-Saint-Gilles", all_data$location), "Haute-Garonne", all_data$location)
all_data$location = ifelse(grepl("Korčula", all_data$location), "Dubrovacko-neretvanska zupanija", all_data$location)
all_data$location = trimws(all_data$location) 

#add NUT3 codes 
codes = match(all_data$location, nuts3_data$location_name)
all_data$NUTS3 = nuts3_data$NUTS3[codes]
all_data = all_data[, c("country", "disease", "year","location", "NUTS3", "case_number", "authors")]
all_data = as.data.frame(all_data)

write.csv(all_data, paste0(path,"Nuts3_cattaneo.csv"), row.names = FALSE)



