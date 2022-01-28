###################################################################################################################
########################### SDOT: Spatial Cell Type Deconvolution by Optimal Transport ############################
## Preliminary work. Under review by the International Conferenceon Machine Learning (ICML). Do not distribute. ###
###################################################################################################################

source("SDOT/helper.R")
source("SDOT/sdot.R")

# Loading sample scRNAseq and spatial data
sample_data <- load_sample_data()


# Estimating the relative abundance of cell types from the sc data
sc_ratios <- table(sample_data$sc_data$labels)/length(sample_data$sc_data$labels)

present_types <- names(sc_ratios)
sc_ratios <- as.numeric(sc_ratios)
names(sc_ratios) <- present_types


# Estimating the signature of cell types from the sc data
sc_signatures <- aggregate(t(as.matrix(sample_data$sc_data$counts)), list(Label = sample_data$sc_data$labels), mean)
rownames(sc_signatures) <- sc_signatures$Label
sc_signatures$Label <- NULL
sc_signatures <- sc_signatures[present_types, ]


# Creating an SDOT object
sdot <- sdot.setup(st_X = t(as.matrix(sample_data$st_data$counts)), 
                   st_coordinates = as.matrix(sample_data$st_data$coordinates), 
                   celltype_signatures = as.matrix(sc_signatures), celltype_ratios = sc_ratios, 
                   st_radius = "auto")

# Solving the model 
sdot.result <- sdot.solve(sdot, iterations = 500)