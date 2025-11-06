#' spin.classifier
#'
#' @param data
#' @param c.logsize.aerosol
#' @param c.depolarization.droplet
#' @param c.logsize.ice
#' @param c.depolarization.ice
#' @param size.ice
#' @param size.all
#' @param processing.cores
#'
#' @return list(data, MODEL.classifier, MODEL.stats, MODEL.pca)
#' @export
#'
##### THIS MODEL USES PARALLEL COMPUTING TO SPEED UP THE ML ALGORITHM

#' This classifier is designed using a three step process
#' 1. Apply preliminary classification using empirically defined values based on control compounds e.g. mineral dust or inorganic sulfate
#' Preliminary classification is based on 20% of the data to initialize the model (two sampling steps)
#' 2. Principle component analysis (PCA) for dimension reduction and identify 95% of the variance
#' 3. Unsupervised support vector machine (SVM) for classification
#' 4. k-Fold Cross Validation
#'
#' @examples
spin.classifier <- function(data,
                            c.logsize.aerosol,
                            c.depolarization.droplet,
                            c.logsize.ice,
                            c.depolarization.ice,
                            size.ice,
                            size.all,
                            processing.cores){

  # Packages needed
  require(dplyr)
  require(doParallel)
  require(caret)
  require(factoextra)

  # Memory intensive so manually force R to clear unused memory
  gc(verbose = F)

  # Set the randomizer seed
  set.seed(123)

  {
    if(missing(c.depolarization.ice)) c.depolarization.ice <- 0.4
    if(missing(c.depolarization.droplet)) c.depolarization.droplet <- 0.2
    if(missing(c.logsize.aerosol)) c.logsize.aerosol <- 0.1
    if(missing(c.logsize.ice)) c.logsize.ice <- 0.4
    if(missing(size.ice)) size.ice <- 2.5
    if(missing(size.all)) size.all <- 0
    if(missing(processing.cores)) processing.cores <- 1
  }

  # -------------------------------------------------------------------------- #
  ##### SECTION: Subset Dataset #####

  {
    # Retrieve indices for thermocouples and bins used for classification
    {
      # Thermocouples
      coldTCs <- colnames(data)[str_which(colnames(data), "C\\d+")]
      warmTCs <- colnames(data)[str_which(colnames(data), "W\\d+")]
    }

    # Create total particle counts for use later
    {
      # Selecting binned data colnames
      bins.ix <- suppressWarnings(which(!is.na(as.numeric(colnames(data)))))
      bins.nm <- colnames(data)[bins.ix]

      # Considering ice as this cutoff size (um) or larger
      bins.ice.nm <- bins.nm[which(as.numeric(bins.nm) >= size.ice)]

      # Total particles
      bins.all.nm <- bins.nm[which(as.numeric(bins.nm) >= size.all)]

      # Calculate total particles above the size limit
      tmp.df <- data %>%
        mutate(`Total Size Ice` = rowSums(select(., .dots = all_of(bins.ice.nm)))) %>%
        mutate(`Total Size All` = rowSums(select(., .dots = all_of(bins.all.nm))))
    }

    # Setup subset dataframe for classification
    tmp.df <- tmp.df %>%
      select(`Lamina Breaks`, `Depolarization`, `OPC Size`, `Lamina S Ice`, `Lamina Temp (C)`, all_of(bins.ix), `Total Size Ice`, `Total Size All`, `Lamina S Liquid`,  all_of(c(coldTCs, warmTCs)), `Inlet Filter ON`, `Warm Wall SP`, `Cold Wall SP`, `Sample Flow (SLPM)`, `Sheath Flow (SLPM)`) %>%
      mutate(`ID` = seq(1:n())) %>%
      mutate(`Class` = "Other") %>%
      mutate(`Lamina Breaks` = as.factor(`Lamina Breaks`)) %>%
      select(`Class`, `ID`, everything())
  }

  # -------------------------------------------------------------------------- #
  ##### SECTION: Preliminary Classification #####

  {
    # AEROSOL
    {
      # This looks at periods where the chamber walls are set for the same temperature and the inlet filter is open
      tmp.df <- tmp.df %>%
        mutate(`Class` = if_else(condition =
                                   (tmp.df$`Warm Wall SP` - tmp.df$`Cold Wall SP` == 0) &
                                   (`OPC Size` <= c.logsize.aerosol) &
                                   (`Depolarization` < c.depolarization.droplet),
                                 true = "Aerosol",
                                 false = `Class`))
    }

    # DROPLET
    {
      tmp.df <- tmp.df %>%
        mutate(`Class` = if_else(condition =
                                   (`OPC Size` >= c.logsize.ice) &
                                   (`Depolarization` < c.depolarization.ice) &
                                   (`Depolarization` >= c.depolarization.droplet) &
                                   (`Total Size Ice` > 0) &
                                   (`Lamina S Liquid` >= 1),
                                 true = "Droplet",
                                 false = `Class`))
    }

    # ICE
    {
      # This looks at periods where the inlet filter is on and 10 um particles are being detected
      # These particles have to be ice from the walls
      tmp.df <- tmp.df %>%
        mutate(`Class` = if_else(condition =
                                   (`OPC Size` >= c.logsize.ice) &
                                   (`Depolarization` >= c.depolarization.ice) &
                                   (`Total Size Ice` > 0),
                                 true = "Ice",
                                 false = `Class`))
    }

    # WATER UPTAKE
    {
      tmp.df <- tmp.df %>%
        mutate(`Class` = if_else(condition =
                                   (`OPC Size` > c.logsize.aerosol) &
                                   (`OPC Size` < c.logsize.ice) &
                                   (`Depolarization` < c.depolarization.droplet) &
                                   (`Lamina S Liquid` < 1),
                                 true = "Water Uptake",
                                 false = `Class`))
    }
  }

  # -------------------------------------------------------------------------- #
  ##### SECTION: Machine Learning #####
  {
    # Setup testing and training datasets
    {
      # Exclude the other category, this is data that is not classified using the linear parameters above
      model.data <- tmp.df %>%
        select(!c(`Warm Wall SP`, `Cold Wall SP`, `Lamina Breaks`, `Lamina S Liquid`, `Total Size Ice`, `Total Size All`)) %>%
        filter(`Class` != "Other") %>%
        mutate(`Class` = factor(`Class`)) %>%
        select(!`ID`)

      # Create a sampling index using caret::createDataPartition
      data.index = caret::createDataPartition(model.data$Class, p = 0.8, list = FALSE)

      # Split data into training and testing
      data.training <- model.data[data.index,]
      data.testing  <- model.data[-data.index,]
    }

    # Parallel machine learning using caret package
    {
      # Repeated K-Fold Cross Validation (10 folds, 3 repeats)
      train_control <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 3)

      # Register the number of cores allowed for processing data
      doParallel::registerDoParallel(cores = processing.cores)

      # Create model classifier
      # Test 15 different parameter combinations
      MODEL.classifier <- caret::train(Class ~ ., data.training, method = "svmRadial", preProcess = c("pca"),
                                       trControl = train_control, tuneLength = 20)

      # Stop parallel cluster
      doParallel::stopImplicitCluster()
    }

    # Retrieve PCA components
    {
      MODEL.pca <- prcomp(data.training[, -1])

      pca.plot <- factoextra::fviz_pca_biplot(MODEL.pca) + factoextra::fviz_screeplot(MODEL.pca, addlabels = T)
    }

    # Test model accuracy
    {
      data.testing$Prediction = predict(MODEL.classifier, data.testing)

      MODEL.stats <- caret::confusionMatrix(data = data.testing$Prediction, reference = as.factor(data.testing$Class))
    }

    # Apply model to main dataset to generate classes
    {
      # Vector of factors
      data$`ML Class` = predict(MODEL.classifier, data)

      data$`Class` = tmp.df$Class

      data <- data %>%
        relocate(c(Class, `ML Class`), .after = `Experiment Ramp ID`)
    }
  }

  # Specify what the output of the function should be
  return.ls <- list(data, MODEL.classifier, MODEL.stats, MODEL.pca)
  return(return.ls)
}
