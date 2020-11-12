################################################################################
#' This function is used to load data into the SOM script
#'
#' @param data_path paths to .nc files. A character vector
#' @param from_date final date
#' @param to_date final date
#' @param var_name names of variables to extract from .nc files. A character vector.
#' @param climatology_name a vector of names of climatology variables. A character vector.
#' @param max_consecutive_na maximum number of consecutive NAs.
#' @param standardize_by if standardization should be by pixel or by date
#' @return a list of dplyr's dataframes
#' @example
#' 
load_data_list <- function(data_path, from_date, to_date, var_name, climatology_name, max_consecutive_na, standardize_by, scaling)
{
    
    # Debug. Following commented lines are used for debug purposes ONLY.
    # rm(list=ls())
    # data_path <- "C:\\Users\\alicemi\\Desktop\\CHR\\SOM_script\\Data"
    # from_date <- "2016-01-01"
    # to_date <- "2016-05-20" # Only for when climatology not required
    # to_date <- "2017-06-01"
    # var_name <- "CHL1_mean"
    # climatology_name <- NA # Only for when climatology not required
    # climatology_name <- "avg_chl"
    # max_consecutive_na <- 2
    # standardize_by <- "pixel"
    # scaling <- "z-score"
    # require(dplyr)
    # require(qchlorophyll)
    
    #-------------------------------------------------------------------------------
    # Loading data
    
    print("Start loading data...")
    
    # Loading data from .nc files
    df <- load_all_as_list(path = data_path,
                           from = from_date,
                           to = to_date,
                           variables = c(var_name)) %>%
        assign_id_and_melt()
    
    #-------------------------------------------------------------------------------
    # Keep coordinates of all the pixels loaded in R
    
    print("Keeping initial reference to coordinates...")
    
    # Keep lon, lat and id_pixel information on loaded data
    df_lon_lat <- df %>% 
        select(all_of(c("lon", "lat", "id_pixel"))) %>%
        distinct()
    
    #-------------------------------------------------------------------------------
    # Remove earth pixels
    
    # Number of pixels used in input
    pixels_IN <- length(unique(df$id_pixel))
    
    v <- as.name(var_name)
    mutate_call <- rlang::exprs(obs_p = sum(!is.na(!!v))/n())
    
    # Remove earth pixels from the analysis
    df <- df %>%
        group_by(id_pixel) %>%
        mutate(!!! mutate_call) %>%
        filter(obs_p > 0) %>%
        ungroup() %>%
        select(-obs_p)
    
    # Percentage of pixels kept for the analysis on the total of pixels used as input
    print(paste("Percentage of pixels kept after earth pixel removal: ", signif(length(unique(df$id_pixel))/pixels_IN * 100, 4), " %", sep = ""))
    
    #-------------------------------------------------------------------------------
    # Stineman interpolation
    
    # Number of rows with all data before interpolation
    complete_rows_before_interpolation <- nrow(df[complete.cases(df), ])
    print(paste("Percentage of rows filled before interpolation: ", signif(complete_rows_before_interpolation/nrow(df) * 100, 4), " %", sep = ""))
    
    # Helper function
    consecutive_NAs <- function(x, n)
    {
        y <- rle(is.na(x))
        y$values <- y$lengths > (n - 0.5) & y$values
        inverse.rle(y)
    }
    
    v_1 <- as.name(var_name)
    mutate_call_1 <- rlang::exprs(var = !! v_1)
    mutate_call_2 <- rlang::exprs(has_NAs_consecutive = consecutive_NAs( !!v_1, !!max_consecutive_na))
    
    # Stineman interpolation
    df <- df %>%
        # Add dummy variable
        mutate(!!! mutate_call_1) %>%
        # Identify groups of consecutive NAs
        group_by(id_pixel) %>%
        mutate(!!! mutate_call_2) %>%
        # Remove grouping
        ungroup() %>%
        
        # Interpolate only pixels which meet requirement
        mutate(var = case_when(has_NAs_consecutive == FALSE ~ imputeTS::na_interpolation(var, option="stine"),
                               TRUE ~ var)) %>%
        # Remove negative pixel from interpolation, if they exist
        #mutate(var = case_when(var < 0 ~ as.numeric(NA),
        #                      TRUE ~ var)) %>%
        # Assign dummy to actual variable
        mutate({{ var_name }} := var) %>%
        # Remove unused variables
        select(-var, -has_NAs_consecutive)
    
    # Number of rows with all data after interpolation
    complete_rows_after_interpolation <- nrow(df[complete.cases(df), ])
    # Calculate percentage of rows filled after interpolation
    print(paste("Percentage of rows filled after interpolation: ", signif(complete_rows_after_interpolation/nrow(df) * 100, 4), " %", sep = ""))
    
    v_2 <- as.name(var_name)
    
    # Perform only if climatology_name is not NA
    if(!is.na(climatology_name))
    {
        #-------------------------------------------------------------------------------
        # Climatology
        
        print("Calculating climatology...")
        
        # Calculate climatology
        df <- df %>%
            # For each pixel, for each id_date
            group_by(id_pixel, id_date) %>%
            # Calculate: climatology, i.e. average value for date, for pixel (avg_chl)
            #           number of observations used in each date to calculate avg_chl (n_observations_used_per_date)
            #           number of dates loaded in R (images per pixel loaded)
            summarise({{climatology_name}} := median({{v_2}}, na.rm = T)) %>%
            # Remove grouping by pixel
            ungroup()
        
        #-------------------------------------------------------------------------------
        # Add back lon, lat to each pixel
        
        print("Adding back lon and lat info...")
        df <- df %>%
            left_join(df_lon_lat, by=c("id_pixel"))
        
    }
    
    #-------------------------------------------------------------------------------
    # From long to wide conversion
    
    print("Conversion from long to wide dataframe...")
    
    if(!is.na(climatology_name))
    {
        # Long to wide conversion
        df <- df %>%
            select(all_of(c("lon", "lat", "id_pixel", "id_date", climatology_name))) %>%
            reshape2::dcast(id_pixel + lat + lon ~ id_date, 
                            value.var = climatology_name)
    }else
    {
        # Long to wide conversion
        df <- df %>%
            select(all_of(c("lon", "lat", "id_pixel", "date", var_name))) %>%
            reshape2::dcast(id_pixel + lat + lon ~ date, 
                            value.var = var_name)
    }
    
    #-------------------------------------------------------------------------------
    # Standardization

    # Select scaling function
    if(scaling == "z-score")
    {
        scaling_function <- scale
        
    }else if(scaling == "min-max")
    {
        scaling_function <- function(x){ (x - min(x, na.rm=T)) / (max(x, na.rm = T) - min(x, na.rm=T)) }
    }else
    {
        print("Normalization function not admitted, using scale...")
        
        scaling_function <- scale
    }
    
    # Standardizing data
    if(standardize_by == "pixel")
    {
        print("Standardizing/scaling by pixel...")
        
        # For each pixel, standardize using mean and sd of the pixel (i.e. standardize by row)
        df_head <- df[, 1:3]
        df_body_names <- names(df[, -c(1:3)])
        df_body <- t(apply(df[, -c(1:3)], 1, scaling_function)) %>% as.data.frame()
        names(df_body) <- df_body_names
        df <- df_head %>%
            bind_cols(df_body)
        
    }else if(standardize_by == "date")
    {
        print("Standardizing/scaling by date...")
        
        # For each date, standardize using mean and sd of the date (i.e. standardize by column)
        df <- df %>%
            as_tibble() %>%
            mutate(across(.cols = !c(id_pixel, lat, lon), .fns = scaling_function))
    }else
    {
        print("Wrong standardization choice... no standardization done...")
    }
    
    # Percentage of pixels kept for the analysis on the total of pixels used as input
    print(paste("Percentage of pixels kept after standardization: ", signif(length(unique(df$id_pixel))/pixels_IN * 100, 4), " %", sep = ""))
    
    #-------------------------------------------------------------------------------
    # Remove NA rows
    df <- df[complete.cases(df), ]
    
    #-------------------------------------------------------------------------------
    # Print number of pixels kept
    
    # Number of pixels kept for the analysis
    pixels_kept <- length(unique(df$id_pixel))
    # Percentage of pixels kept for the analysis on the total of pixels used as input
    print(paste("Percentage of pixels without any NA kept: ", signif(pixels_kept/pixels_IN * 100, 4), " %", sep = ""))
    
    #-------------------------------------------------------------------------------
    # Output
    
    # Output to be returned
    list_df <- list(wide_data = df, original_grid_df = df_lon_lat)
    
    # Return
    return(list_df)
}
