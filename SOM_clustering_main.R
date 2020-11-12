# Script version

# V. 0.6

#-------------------------------------------------------------------------------
# This script performs the following:

# 1. Load of data
# 2. Interpolation of missing data
# 3. Conversion of data format from long to wide
# 2. Calculates best SOM size
# 3. Gets prototypes from best SOM and clusters for each pixel through the prototypes
# 4. Clustering of prototypes (hierarchical clustering)
# 5. Associate to each pixel the cluster of the prototypes:
#       pixel cluster belongs to a prototype (A) ->
#       prototype is clustered (B)
#   From (B) to (A)
#   

#-------------------------------------------------------------------------------
# Clear environment

rm(list = ls())

#-------------------------------------------------------------------------------
# Load required libraries

require(qchlorophyll)   # V. >= 2.1
require(dplyr)          # V. >= 1.0.2
require(SOMbrero)       # V. >= 1.3-1
#require(fdm2id)        # V. >= 0.9.3
#require(reshape2)      # V. >= 1.4.4
require(mclust)         # V. >= 5.4.6
require(factoextra)     # V. >= 1.0.7
#require(readr)         # V. >= 1.3.1
#require(imputeTS)      # V. >= 3.1
#require(rlang)         # V. >= 0.4.7

#-------------------------------------------------------------------------------
# Load auxiliary functions

# Auxiliary functions path
AUX_FUNCTIONS_PATH <- "C:\\Users\\alicemi\\Desktop\\CHR\\SOM_script\\Source\\auxiliary_functions"

# Load functions
source(paste(AUX_FUNCTIONS_PATH, "load_data_list.R", sep = "\\"))

# Clean workspace
rm(AUX_FUNCTIONS_PATH)
#-------------------------------------------------------------------------------
# Set script parameters

# Initial date of nc files
FROM_DATE <- "2016-01-01"
# Final date of nc files
TO_DATE <- "2017-01-01"

# Note:
# The length of:
# - DATA_PATH
# - VAR_NAME
# - MAX_CONSECUTIVE_NA
# - CLIMATOLOGY_NAME (only if climatology is required)
#
# must always be the same. CLIMATOLOGY_NAME can also be NA in case climatology is not required. 

# Paths to nc files
DATA_PATH <- c("C:\\Users\\alicemi\\Desktop\\CHR\\SOM_script\\Data",
               "C:\\Users\\alicemi\\Desktop\\CHR\\SOM_script\\Data")

# Variables to extract from nc files
VAR_NAME <- c("CHL1_mean",
              "CHL1_mean")

# Names to assign to climatology. Put NA if climatology is not required
CLIMATOLOGY_NAME <- c("avg_chl",
                      "avg_chl_2")

# Maximum number of consecutive NAs to be allowed (example: if it is equal to 3, only gaps of 2 or less NA in the pixel time series will be interpolated)
MAX_CONSECUTIVE_NA <- c(5,
                        5)

# Som size vector (only square matrices allowed: this will try 2x2, ..., 5x5)
SOM_SIZE_VECTOR <- 2:5

# Standardization criterion (possible values are "pixel" or "date" to standardize by pixel or by date)
STANDARDIZE_BY <- "pixel"

# Scaling function to be used. Other acceptable values are: min-max
SCALING <- "z-score"

# Higher and lower percentile of topographic error for SOM selection
LOWER_PERCENTILE_SOM_SELECTION <- 0.0
HIGHER_PERCENTILE_SOM_SELECTION <- 0.4

# Path of Log
LOG_PATH <- "C:\\users\\alicemi\\Desktop\\CHR\\SOM_script"

# Log name
LOG_NAME <- "log.csv"

# Criterion for selecting best cluster. Cluster closest to 90th percentile interclass inertia value is the best. This value should be between 0 and 1
CLUSTER_PERCENTILE_SELECTION <- 0.90

# Path where to save the results
OUTPUT_PATH <- "C:\\users\\alicemi\\Desktop\\CHR\\SOM_script"

# Output dataframe name
OUTPUT_DF_NAME <- "SOM_CLUSTERING_CLASSIFICATION.csv"

# Set seed for reproducibility purposes
SEED <- 10

#-------------------------------------------------------------------------------
# Load required variables from .nc files

# List containing loaded data
data_list <- list()

# Start loading data for each set of paths
for(i in 1:length(DATA_PATH))
{
    # Print progress info
    print(paste("Loading data variable: ", VAR_NAME[i]))
    # Fill list
    data_list[[i]] <- load_data_list(DATA_PATH[i],
                                     FROM_DATE,
                                     TO_DATE,
                                     VAR_NAME[i],
                                     CLIMATOLOGY_NAME[i],
                                     MAX_CONSECUTIVE_NA[i],
                                     STANDARDIZE_BY,
                                     SCALING)
}

rm(DATA_PATH, FROM_DATE, TO_DATE, MAX_CONSECUTIVE_NA, STANDARDIZE_BY, i, load_data_list, SCALING)
#-------------------------------------------------------------------------------
# Keep the initial dataframe of all the pixels loaded in R

# initial_pixel_df contains all the pixels loaded in R
initial_pixel_df <- data_list[[1]]$original_grid_df
data_list[[1]]$original_grid_df <- NULL

# All the pixels initially loaded in R are contained in initial_pixel_df dataframe
# All the pixels that will be analyzed by the SOM algorithm are contained in min_common_den dataframe

#-------------------------------------------------------------------------------
# Perform the analysis only on the common base of pixels so that SOM algorithm can work without NAs

# Common base of pixels used in the analysis
min_common_den <- data.frame(id_pixel = data_list[[1]]$wide_data$id_pixel) %>%
    as_tibble()

# If there is more than 1 element in data list
if(length(data_list) > 1)
{
    # For each element in data list
    for(i in 2:length(data_list))
    {
        # Extract id_pixel
        temp <- data.frame(id_pixel = data_list[[i]]$wide_data$id_pixel) %>%
            as_tibble()
        # Make inner join to keep only common pixels
        min_common_den <- min_common_den %>%
            inner_join(temp, by="id_pixel")
        # Delete duplicated df
        data_list[[i]]$original_grid_df <- NULL
    }
    
    # Now min_common_den contains only common pixels to all the loaded data
    
    # Remove temporary df
    rm(temp)
}

# For each wide dataframe, keep only common pixels
for(i in 1:length(data_list))
{
    data_list[[i]]$wide_data <- data_list[[i]]$wide_data[data_list[[i]]$wide_data$id_pixel %in% min_common_den$id_pixel, ]
}

rm(i)
#-------------------------------------------------------------------------------
# Bind by column the data to be fed to the SOM algorithm

# Assign unique names to columns of som_df. Account for when climatology is NA
if(!is.na(CLIMATOLOGY_NAME[1]))
{
    col_names_vector <- CLIMATOLOGY_NAME
}else
{
    col_names_vector <- VAR_NAME
}

# Dataframe to be fed to the SOM algorithm
som_df <- data_list[[1]]$wide_data[2:ncol(data_list[[1]]$wide_data)]
# Assign unique names to the columns
names(som_df) <- c(names(som_df)[1:2], paste(col_names_vector[1], names(som_df[3:ncol(som_df)]), sep = "_"))

# Perform the same operation on the other variables of data list (if there is more than 1)
if(length(data_list) > 1)
{
    # For each variable in data list
    for(i in 2:length(data_list))
    {
        # Assign data to a temporary dataframe
        temp <- data_list[[i]]$wide_data[, 4:ncol(data_list[[i]]$wide_data)]
        # Assign unique name to the columns
        names(temp) <- paste(col_names_vector[i], names(temp), sep = "_")
        # Bind by column the data to som_df
        som_df <- cbind(som_df, temp)
    }
    # Remove temporary dataframe
    rm(temp)
}

rm(CLIMATOLOGY_NAME, VAR_NAME, i, data_list, col_names_vector)
#-------------------------------------------------------------------------------
# SOM algorithm

# Set seed for reproducibility purposes
set.seed(SEED)

# List of topographic error
topo_error <- list()
# List of quantization error
quant_error <- list()
# List of calculated SOM
SOM_list <- list()

# Start to calculate SOM output for each one of the selected sizes
k <- 1
for(i in SOM_SIZE_VECTOR)
{
    # Print info on progress
    print(paste("Training SOM of size: ", i, "x", i, "...", sep=""))
    
    # Calculate SOM
    som_out <- trainSOM(som_df, dimension = c(i, i),
                        topo = "hexagonal",
                        radius.type = "gaussian",
                        dist.type = "euclidean",
                        maxit = 5000,
                        init.proto = "obs",
                        scaling = "none",
                        nb.save = 25)
    
    # Calculate SOM quality
    qsom <- quality(som_out)
    
    # Print info on SOM quality
    print(paste("Topographic error: ", qsom[1], " || Quantization error: ", qsom[2]))
    
    # Assign calculated quantities to the respective list
    topo_error[[k]] <- qsom$topographic
    quant_error[[k]] <- qsom$quantization
    SOM_list[[k]] <- som_out
    
    k <- k + 1
}

# Conversion of quantization error from list to vector
quant_error <- unlist(quant_error)
# Conversion of topographic error from list to vector
topo_error <- unlist(topo_error)

rm(i, k, qsom, som_out, SEED)
#-------------------------------------------------------------------------------
# Best SOM size selection

# Criterion:
# 1. Calculate QE and TE for each SOM
# 2. Select SOM with TE between 5th and 10th percentiles
# 3. Within these, select SOM with lowest QE

# Quantiles of topographic error
quantiles <- quantile(topo_error, c(LOWER_PERCENTILE_SOM_SELECTION, HIGHER_PERCENTILE_SOM_SELECTION))

# Vector of SOM in quantiles
som_in_quantiles <- between(topo_error, quantiles[1], quantiles[2])

if(sum(som_in_quantiles) == 1)
{
    # Print only 1 SOM size matches the selected percentile interval
    som_best_size <- SOM_SIZE_VECTOR[som_in_quantiles]
    
    
}else if(sum(som_in_quantiles) > 1)
{
    
    # More than 1 SOM size matches the selected percentile interval
    som_best_size <- SOM_SIZE_VECTOR[som_in_quantiles][which.min(quant_error[som_in_quantiles])]
    
}else
{
    # No SOM size matches the selected percentile interval
    print("No SOM matches the percentiles given... Try again increasing percentiles..")
}

# Print info
print(paste("SOM size which minimizes the total error: ", som_best_size, "x", som_best_size, sep=""))

# Get SOM with best quality
best_som_out <- SOM_list[[which(SOM_SIZE_VECTOR == som_best_size)]]

# Plot energy
plot(best_som_out, what = "energy")

rm(SOM_list, quantiles, som_in_quantiles, LOWER_PERCENTILE_SOM_SELECTION, HIGHER_PERCENTILE_SOM_SELECTION)
#-------------------------------------------------------------------------------
# Save log

print(paste("Saving log to: ", LOG_PATH, "\\", LOG_NAME, sep = ""))

# Dataframe to log
log_df <- data.frame(SOM_size_squared = SOM_SIZE_VECTOR,
                     topo_error = topo_error,
                     quantization_error = quant_error,
                     SOM_size_selection = som_best_size)

# Write dataframe to .csv file
readr::write_csv(log_df, path = paste(LOG_PATH, LOG_NAME, sep="\\"))

rm(topo_error, quant_error, log_df, LOG_NAME, LOG_PATH, SOM_SIZE_VECTOR)
#-------------------------------------------------------------------------------
# Select cluster size

# Acceptable cluster sizes go from 2 to SOM size^2
cluster_size_vector <- 2:som_best_size^2

rm(som_best_size)
#-------------------------------------------------------------------------------
# Clustering of the prototypes

# List of prototypes cluster quality
cluster_quality <- list()

k <- 1
for(j in cluster_size_vector)
{
    # Hierarchical clustering of the prototypes
    clustering_out <- superClass(best_som_out,
                                 method = "complete",
                                 k = j)
    
    # Clusters
    cl <- clustering_out$cluster
    
    # Calculate cluster quality
    cluster_quality[[k]] <- fdm2id::intern.interclass(cl, clustering_out$som$prototypes)
    
    k <- k + 1
}

# Convert cluster_quality from list to vector
cluster_quality <- unlist(cluster_quality)

rm(j, k, cl)
#-------------------------------------------------------------------------------
# Choose best cluster

# Plot hitmap
plot(clustering_out, what="obs", type="hitmap")

plot(cluster_size_vector, cluster_quality, type="l", xlab="Number of clusters k", ylab="Interclass inertia", main = "Interclass inertia", lwd = 2)
points(cluster_size_vector, cluster_quality)

# Line of 90% of max value
abline(a = max(cluster_quality)*CLUSTER_PERCENTILE_SELECTION, b=0, col="blue", lwd = 3)

# Change in interclass inertia
plot(cluster_size_vector[1:(length(cluster_quality) - 1)], diff(cluster_quality), col = "red", lwd = 2, type = "l", xlab = "Cluster size", ylab = "Change in interclass inertia", main = "Change in interclass inertia")
points(cluster_size_vector[1:(length(cluster_quality) - 1)], diff(cluster_quality), col = "red")

# Select best cluster size. Best cluster is the one closest to the 90th percentile value of interclass cluster inertia
best_cluster_size <- cluster_size_vector[which(min(abs(cluster_quality - max(cluster_quality)*CLUSTER_PERCENTILE_SELECTION)) == abs(cluster_quality - max(cluster_quality)*CLUSTER_PERCENTILE_SELECTION))]

rm(CLUSTER_PERCENTILE_SELECTION)
#-------------------------------------------------------------------------------
# Get best clustering result

# Get best clustering result
clustering_out <- superClass(best_som_out,
                             method = "complete",
                             k = best_cluster_size)
# Perform Anova
summary(clustering_out)

rm(cluster_quality, cluster_size_vector, best_cluster_size)
#-------------------------------------------------------------------------------
# Validation using Mclust package

# Perform clustering
mc <- Mclust(best_som_out$prototypes)

# Summary of the clustering output
summary(mc)
# Optimal selected model name
mc$modelName
# Optimal number of clusters
mc$G
# Probability to belong to a given cluster
head(mc$z, 5)

#-------------------------------------------------------------------------------
# Plot the results

# BIC values used for choosing the number of clusters
fviz_mclust(mc, "BIC", palette = "jco")
# Classification: plot showing the clustering
fviz_mclust(mc, "classification", geom = "point", 
            pointsize = 1.5, palette = "jco")
# Classification uncertainty
fviz_mclust(mc, "uncertainty", palette = "jco")

#-------------------------------------------------------------------------------
# Getting data back to original grid

# SOM assigns to each pixel a cluster by assigning it to a prototype (neuron)
# Hierarchical clustering forms clusters of prototypes
# Now we need to assign to each pixel the cluster of the corresponding prototype

# For example:
# With a 5x5 SOM:
#   SOM algorithm assigns each pixel to one of 5x5=25 prototypes (neurons)
#   Hierarchical clustering assigns a cluster to each prototype for example, with k=3 we have 3 total cluster
#   Now each prototype will belong to one of these 3 clusters
#   Each pixel must now be assigned to the same cluster of the prototype it belongs to.

# Assign to each pixel analyzed by the SOM the cluster (prototype) to which it is assigned by the SOM
som_df <- som_df %>%
    mutate(prot_id = best_som_out$clustering)

# Associate to each prototype an id and the cluster calculated via hierarchical clustering.
proto_df <- data.frame(prot_id = 1:nrow(clustering_out$som$prototypes)) %>%
    as_tibble() %>%
            # Result from SOMbrero package
    mutate(cluster_proto = clustering_out$cluster,
           # Result from Mclust package
           cluster_proto_mclust = mc$classification)

# Associate to each pixel the cluster of the prototype to which they belong
som_df <- som_df %>%
    left_join(proto_df, by="prot_id")

rm(mc, best_som_out)
#-------------------------------------------------------------------------------
# Get back clustering result to initial dataframe

# Associate to pixel in the analysis the final cluster
min_common_den <- min_common_den %>%
            # Result from SOMbrero package
    mutate(cluster = som_df$cluster_proto,
           # Result from Mclust package
           cluster_mclust = som_df$cluster_proto_mclust)

# Copy the final cluster in the initial dataframe loaded in R. This is the output of the script
initial_pixel_df <- initial_pixel_df %>%
    left_join(min_common_den, by="id_pixel")

rm(clustering_out, min_common_den, proto_df, som_df)
#-------------------------------------------------------------------------------
# Save output

# Save initial pixel df
readr::write_csv(initial_pixel_df, path = paste(OUTPUT_PATH, OUTPUT_DF_NAME, sep="\\"))

# Print debug info
print(paste("Output dataframe saved to: ", OUTPUT_PATH, OUTPUT_DF_NAME, sep=""))

rm(OUTPUT_DF_NAME, OUTPUT_PATH)
#-------------------------------------------------------------------------------
# End of script

print("Done!")
