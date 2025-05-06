##### SCRIPT: SPectrometer for Ice Nucleation (SPIN) Level 1 #####
#' @author Christopher Rapp
#' @description
#'
#'
#'
# ---------------------------------------------------------------------------- #
##### SECTION: Libraries, Functions, Paths #####
#' These are directories that are in use for this code.
#' If files or scripts are renamed, moved, or changed, the code will NOT work!
#' Please keep updated for future readers of this code!
#'

{
  # Clears global environment and all saved variables
  rm(list = ls())
  gc(verbose = F)

  # Data wrangling libraries
  library(data.table)
  library(lubridate)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(collapse)
  library(purrr)

  # Plotting libraries
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(ggsci)
  library(viridis)
  library(RColorBrewer)

  # ML libraries
  library(caret)
  library(factoextra)

  # Plotting parameters
  resolution.dpi = 400
  font.family = "Helvetica"

  # Set working directory
  setwd("~/Library/CloudStorage/Box-Box/BVOC Chamber Study/")
  work.dir <- getwd()

  # Instrument time zone
  tz.c = "US/Eastern"

  #' @import
  #' Specify import directory
  import.spin = paste0(work.dir, "/SPIN/export/level0/")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/SPIN/export/level1/")
  export.plot = paste0(work.dir, "/SPIN/export/plots/")

  #' @importFrom
  #' Function scripts
  source(paste0("~/Documents/GitHub/functions/", "functions_microphysics.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_plotting.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_lamina.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_classification.R"))
}

# ---------------------------------------------------------------------------- #
##### SECTION: File Selection #####
#'

{
  # Directories
  spin.dirs <- list.dirs(path = import.spin,
                         recursive = FALSE,
                         full.names = FALSE)

  # Pull POSIX style dates from directory names
  spin.dates <- as.Date(spin.dirs, format = '%Y%m%d')

  # Generate file paths
  spin.path  <- paste0(import.spin, spin.dirs)
}

# ---------------------------------------------------------------------------- #
##### SECTION: Level 1 #####
#'

# Plotting on?
PLOT.ON = T

{
  {
    # Set threshold for size total
    THRESHOLD.Dp.nm <- 2.5
    data.export.ls <- NULL
    model.export.ls <- NULL
    for (n in 1:length(spin.dates)){

      # Code Timing
      pct <- proc.time()

      {
        # File selection and pre-processing
        {
          print(paste0("Level 1 SPIN Data for ", spin.dates[n], " from directory ", spin.path[n]))

          # Memory management
          gc(verbose = F)

          if (spin.dates[n] == "2024-01-24"){
            next
          }

          # List all SPIN files exported
          files.spin <- list.files(
            path = spin.path[n],
            recursive = TRUE,
            full.names = TRUE,
            pattern = '*.csv'
          )

          if (length(files.spin) != 3){
            next
          }
        }

        # -------------------------------------------------------------------- #
        ##### SUBSECTION: SPIN Data #####
        #' Data for SPIN is split into three files
        #' The first is the main file that contains housekeeping data and binned counts
        #' The second is the log file which is used to determine when certain sequences have started
        #' The third is the particle by particle data recorded by the OPC
        #'

        {
          {
           # Loop through SPIN files and read in
            data.ls <- lapply(files.spin, function(x){
              tmp.df <- fread(x)
              return(tmp.df)
            })

            # Retrieve list components
            {
              dataLOG.df <- data.ls[[which(lengths(data.ls) <= 3)]]
              rawPBP.df <- data.ls[[which(lengths(data.ls) == 8)]]
              rawSPIN.df <- data.ls[[which(lengths(data.ls) > 100)]]
            }

            # Formatting time strings
            dataLOG.df <- dataLOG.df %>%
              mutate(`Local Time` = lubridate::with_tz(`Local Time`, tzone = tz.c))

            # Retrieve experiment start time from LOG file
            sequence.ramps.tm = dataLOG.df$`UTC Time`[str_which(dataLOG.df$Message, "Starting Sequence Update SPs from Ramp.")]

            # Skip loops with no ramp (indicating experiment never started)
            {
              ramp.break <<- FALSE

              tryCatch(if (length(sequence.ramps.tm) < 1) {stop("Error")},
                       error = function(e) {ramp.break <<- TRUE})

              # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
              if (ramp.break) {next}
            }

            # Skip loops with empty PBP data
            {
              PbP.break <<- FALSE

              tryCatch(if (nrow(rawPBP.df) == 0) {stop("Error")},
                       error = function(e) {PbP.break <<- TRUE})

              # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
              if (PbP.break) {next}
            }

            # Filter data after an experimental ramp has started
            dataSPIN.df <- rawSPIN.df %>%
              filter(`UTC Time` >= sequence.ramps.tm[1]) %>%
              mutate(`Local Time` = lubridate::with_tz(`Local Time`, tzone = tz.c)) %>%
              mutate(`Start Time` = lubridate::with_tz(`Local Time`, tzone = tz.c))

            # The timestamp variable is sometimes converted to int64
            # Reason is unknown
            dataSPIN.df$Timestamp <- as.numeric(dataSPIN.df$Timestamp)
          }

          # ------------------------------------------------------------------ #
          ##### SUBSECTION: Binned Data #####

          # Bin formatting
          {
            bins.ix <- suppressWarnings(which(!is.na(as.numeric(colnames(dataSPIN.df)))))
            bins.nm <- colnames(dataSPIN.df)[bins.ix]

            # Subset bins
            # Drop the first as it is not defined by a lower limit
            dataBINS.df <- dataSPIN.df %>%
              select(all_of(bins.ix)) %>%
              setnames(., new = paste0(bins.nm))

            # ---------------------------------------------------------------- #
            ##### SUBSECTION: Background Removal #####
            {
              # Considering ice as this cutoff size (um) or larger
              cutoff.bins.nm <- bins.nm[which(as.numeric(bins.nm) >= THRESHOLD.Dp.nm)]

              # Calculate total particle density larger than cutoff
              tmp.df <- dataBINS.df %>%
                mutate(`Total Size Ice` = rowSums(select(., .dots = all_of(cutoff.bins.nm)))) %>%
                mutate(`Total Size All` = rowSums(select(., .dots = all_of(bins.nm)))) %>%
                mutate(`Inlet Filter ON` = dataSPIN.df$`Inlet Filter ON`) %>%
                select(`Inlet Filter ON`, `Total Size Ice`, `Total Size All`)

              # Set total backgrounds to NA for when the inlet filter is on
              # Using 0 will sway statistics significantly
              tmp.df <- tmp.df %>%
                mutate(`Total Background` = if_else(`Inlet Filter ON` == 1, `Total Size Ice`, NA), .after = `Total Size Ice`)

              # Using a linear interpolation to fill in the time series than averaging gives a higher background ice concentration
              # This gives a more conservative approach in declaring ice counts
              tmp.df <- tmp.df %>%
                mutate(`Total Background` = zoo::na.approx(`Total Background`, na.rm = F))

              # Use a rolling average
              tmp.df <- tmp.df %>%
                mutate(`Total Background` = ceiling(zoo::rollapplyr(`Total Background`, 5*60, mean, fill = 0)))

              # Round up to the nearest particle
              tmp.df <- tmp.df %>%
                mutate(`Total Background` = ceiling(`Total Background`)) %>%
                select(!`Inlet Filter ON`) %>%
                mutate(`Size Threshold` = THRESHOLD.Dp.nm)

              # Merge totals back including backgrounds
              dataBINS.df <- cbind(dataBINS.df, tmp.df)
            }

            # Conversion from counts to number density
            {
              # Unit conversion from counts to concentration
              # Convert flow rate to cm^3 of air per second
              # Use the sample flow not the sheath as the OPC is designed to only look at the focused lamina
              convert.nm = (dataSPIN.df$`Sample Volumetric Flow (LPM)`*1000)/60

              # Axis 1 indicates perform rowwise
              dataBINS.df <- sweep(dataBINS.df, 1, convert.nm, "/")
            }

            # Merge this back to the main dataframe replacing the counts
            dataSPIN.df <- dataSPIN.df %>%
              select(!all_of(bins.ix)) %>%
              cbind(dataBINS.df) %>%
              mutate(`Units` = "n/cc") %>%
              relocate(all_of(bins.nm), .before = "Units") %>%
              relocate(c(`Total Size Ice`, `Total Background`), .after = `Units`)
          }

          # ------------------------------------------------------------------ #
          ##### SUBSECTION: PBP Data #####

          {
            # Low and high gain amplifier offsets
            rawPBP.df$S1[rawPBP.df$S1 > 14930] = rawPBP.df$S1[rawPBP.df$S1 > 14930] - 2000
            rawPBP.df$P1[rawPBP.df$P1 > 13570] = rawPBP.df$P1[rawPBP.df$P1 > 13570] - 2000
            rawPBP.df$P2[rawPBP.df$P2 > 14610] = rawPBP.df$P2[rawPBP.df$P2 > 14610] - 3000

            # Benchmarking
            x1 <- nrow(rawPBP.df)

            # Very very very fast method to do averaging for this dataset
            rawPBP.df <- rawPBP.df %>%
              collapse::fgroup_by(`Time Stamp`) %>%
              collapse::fmean()

            # Benchmarking
            x2 <- nrow(rawPBP.df)

            # Remove problematic values
            rawPBP.df <- na.omit(rawPBP.df)
            rawPBP.df <- rawPBP.df[!is.infinite(rowSums(rawPBP.df)), ]

            # Benchmarking
            print(paste0(spin.dates[n], ": ", x1, " PbP observations averaged to ", x2))
            rm(x1, x2)

            # Calculate log values
            rawPBP.df$`Log P1` <- log10(rawPBP.df$P1)
            rawPBP.df$`Log P2` <- log10(rawPBP.df$P2)
            rawPBP.df$`Log S1` <- log10(rawPBP.df$S1)
            rawPBP.df$`Log Size` <- log10(rawPBP.df$Size)

            # Remove any problematic values caused by taking the log
            rawPBP.df <- rawPBP.df[!is.nan(rowSums(rawPBP.df)), ]
            rawPBP.df <- rawPBP.df[!is.infinite(rowSums(rawPBP.df)), ]
            rawPBP.df <- na.omit(rawPBP.df)

            # Subset
            dataPBP.df <- rawPBP.df %>%
              select(`Time Stamp`, `S1`, `P1`, `P2`, `Size`, `Log S1`, `Log P1`, `Log P2`, `Log Size`)

            # Change name for merging later
            setnames(dataPBP.df, old = "Time Stamp", new = "Time (sec)")
          }

          # ------------------------------------------------------------------ #
          ##### SUBSECTION: Aggregate Data #####

          {
            # Time stamps need to be averaged to the nearest hundreth second to merge properly
            dataSPIN.df$`Time (sec)` <- round(dataSPIN.df$`Time (sec)`, 1)
            dataPBP.df$`Time (sec)` <- round(dataPBP.df$`Time (sec)`, 1)

            # Perform joins
            # Right join is focused on only allowing observations from the first join the second if there is a match
            # Inner join is most conservative and only allows observations that match in both
            dataALL.df <- left_join(dataSPIN.df, dataPBP.df, by = "Time (sec)")

            # Filter out any empty time values that occur when merging datasets
            dataALL.df <- dataALL.df %>%
              filter(!is.na(`UTC Time`))

            # Date string for plotting later
            date.c <- unique(lubridate::as_date(dataALL.df$`UTC Time`))[1]

            # Find lamina breaks that segment dataset
            dataALL.df <- dataALL.df %>%
              mutate(`Lamina Breaks` = round(`Lamina Temp (C)`/50, 1)*50)

            # Calculate the depolarization ratio
            # Normalize values to be between 0 and 1 for standardization
            dataALL.df <- dataALL.df %>%
              mutate(`Pseudo Depolarization` = `S1`/(`P1` + `P2`)) %>%
              mutate(`Pseudo Depolarization` = if_else(`Pseudo Depolarization` > 1, NA, `Pseudo Depolarization`)) %>%
              mutate(`Pseudo Depolarization` = if_else(`Pseudo Depolarization` < 0, NA, `Pseudo Depolarization`)) %>%
              mutate(`Depolarization` = 2*`Pseudo Depolarization`) %>%
              mutate(`Depolarization` = scales::rescale(`Depolarization`))

            # Scale the Log Size variable so it is between 0 and 1
            dataALL.df <- dataALL.df %>%
              mutate(`OPC Size` = scales::rescale(`Log Size`, to = c(0, 1)), .before = `Depolarization`)

            # Remove NA values for depolarization as these particles should not be considered for analysis
            dataALL.df <- dataALL.df %>%
              filter(!is.na(`Depolarization`)) %>%
              filter(!is.na(`Pseudo Depolarization`))

            # Scaling necessary to convert from the raw depolarization uisng the sum to the proper equation with the mean
            {
              # Find the slope to match the new normalized data back to the original data
              scaling.factor.ml <- lm(formula = `Pseudo Depolarization` ~ `Depolarization`, data = dataALL.df)

              # Applying the linear regression correction factor
              # This is simply a slope adjustment
              dataALL.df <- dataALL.df %>%
                mutate(`Depolarization` = `Depolarization`*coefficients(scaling.factor.ml)[2])

              # Plot path
              export.plot.path = paste0(export.plot, str_remove_all(date.c, '-'), "/")

              plot.filename <- paste0(export.plot.path, date.c, " Depolarization Scaling.png")

              corr.R = signif(cor(dataALL.df$`Pseudo Depolarization`, dataALL.df$`Depolarization`, use = "complete.obs"), 6)

              label = paste0("R^2*' = '*", corr.R)

              correlation.check.gg <- ggplot(dataALL.df, aes(x = `Pseudo Depolarization`, y = Depolarization)) +
                geom_point(aes(x = `Pseudo Depolarization`, y = `Depolarization`/coefficients(scaling.factor.ml)[2], color = "Normalized")) +
                geom_point(aes(x = `Pseudo Depolarization`, y = `Depolarization`, color = "Normalized/Scaled")) +
                geom_abline(slope = 1) +
                annotate("text", x = 0.05, y = 0.95, label = label, parse = T) +
                labs(title = date.c, subtitle = "Scaling Conversion") +
                theme_minimal() +
                coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "on") +
                theme(legend.title = element_blank())

              ggsave(filename = plot.filename, correlation.check.gg, bg = "white", width = 7, height = 6)
            }

            # Order columns
            dataALL.df <- dataALL.df %>%
              relocate(`S1`:`Depolarization`, .after = `Inlet Filter ON`)

            print(paste0(date.c, ": Data Aggregation Complete"))

            # Clean up memory
            rm(rawPBP.df, rawSPIN.df, correlation.check.gg, scaling.factor.ml, corr.R)
            gc(verbose = F)
          }

          # ------------------------------------------------------------------ #
          ##### SUBSECTION: Filtering Data #####

          {
            # Standard error removal of irregular flows
            # This can happen during an experiment where flow upstream was plugged
            # Only applied if there are multiple large shifts in flow to prevent unncessarily applying a filter to good data
            flows.nm <- abs(diff(dataALL.df$`Total Flow (LPM)`))

            if (length(which(flows.nm > 1)) > 3){

              # Outlier detection
              threshold.lower.nm = quantile(dataALL.df$`Total Flow (LPM)`, probs = c(0.25)) - 1.5*IQR(dataALL.df$`Total Flow (LPM)`)
              threshold.upper.nm = quantile(dataALL.df$`Total Flow (LPM)`, probs = c(0.75)) + 1.5*IQR(dataALL.df$`Total Flow (LPM)`)

              # Remove data outside the outliers
              dataALL.df <- dataALL.df %>%
                filter(`Total Flow (LPM)` >= threshold.lower.nm) %>%
                filter(`Total Flow (LPM)` <= threshold.upper.nm)
            }
          }

          # ------------------------------------------------------------------ #
          ##### SUBSECTION: Classification #####
          # This needs to be performed on a experiment by experiment basis
          # Days were two experiments were conducted need to be split for classification

          {
            print(paste0(date.c, ": Classification"))

            # Split dataframe using base R functionality
            tmp.ls <- split(dataALL.df, dataALL.df$`Experiment Ramp ID`)

            # Loop through split dataframe and run classifier for each
            data.ls <- NULL
            index.ls <- NULL
            for (i in 1:length(tmp.ls)){

              # Classifier with error control
              {
                SVM.break <<- FALSE

                # Try applying the classifier
                tryCatch(expr = {

                  # Apply classifier function
                  # See documentation for specifics
                  data.ls[[i]] <- spin.classifier(data = tmp.ls[[i]],
                                                  c.logsize.aerosol = 0.125,
                                                  c.logsize.ice = 0.4,
                                                  c.depolarization.droplet = 0.16,
                                                  c.depolarization.ice = 0.4,
                                                  size.ice = 2.5,
                                                  size.all = 0,
                                                  processing.cores = 12)

                  # Index keeping for experiments and dates
                  data.ls[[i]] <- append(data.ls[[i]], list(date.c, unique(tmp.ls[[i]]$`Experiment Ramp ID`), "PASS"))

                }, error = function(e){
                  SVM.break <<- TRUE
                })

                # Stop code from continuing
                if (SVM.break) {

                  # Index keeping for experiments and dates
                  data.ls[[i]] <- append(c(NA, NA, NA, NA), list(date.c, unique(tmp.ls[[i]]$`Experiment Ramp ID`), "FAIL"))

                  next
                }
              }
            }

            # Remove any iterations at which there is a NA rather than a model run
            tmp.ls <- compact(map(data.ls, 1))

            if (all(is.na(map(tmp.ls, 1)))){

              next
            }

            # Remove empty lists for merging if a failed experiment occured
            tmp.ls <- tmp.ls[!is.na(tmp.ls)]

            # Merge original dataframe back from classifier
            export.df <- rbindlist(tmp.ls)

            # Merge all model type data to a list while retaining structure
            model.ls <- lapply(data.ls, "[", 2:7)

          }
        }
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Export Data #####
      #'
      #'

      if (nrow(export.df) != 0){

        # Export data to a temporary list for optional plotting later
        data.export.ls[[n]] <- export.df

        # Export modeling data to an optional list
        model.export.ls[[n]] <- model.ls

        # Check if export path exists
        # If it does not, create it
        if (!dir.exists(export.data)) {
          # Create a dated directory to send plots to
          dir.create(export.data, mode = "777")
        }

        export.filename = paste0(export.data, "SPIN001_SPIN_", str_remove_all(date.c, '-'), "_level1.csv")

        print(paste0("Exporting Level 1 Data: ", export.filename))

        # Save data using data.table::fwrite
        data.table::fwrite(export.df, file = export.filename, showProgress = T)

        # Code benchmarking
        print(paste0(date.c, ', Elapsed Time (s) ', round((
          proc.time() - pct)[3], 2)))
      }
    } # END LOOP

    # Remove empty data runs
    data.export.ls <- compact(data.export.ls)
    model.export.ls <- compact(model.export.ls)

    # Backup data from current run for testing
    saveRDS(data.export.ls, file = paste0(export.data, "TMP.RDS"))

    # Backup modeling data
    saveRDS(model.export.ls, file = paste0(export.data, "MOD.RDS"))

    # Clean environment
    rm(data.ls, dataALL.df, dataLOG.df, dataPBP.df, dataSPIN.df, bins.ix, bins.nm, convert.nm, date.c, n, ramp.break, sequence.ramps.tm)
    gc(verbose = F)
  }

  # -------------------------------------------------------------------------- #
  ##### SECTION: Plotting #####
  #' Classification
  #' Background Removal
  #' Correlation Plots
  #' Probability Density of Univariate PBP Data
  #' Bivariate KDE Joint Probability Density
  #'

  if (PLOT.ON == T){

    # # Read back in looped data from above and filter empty data
    {
      # Dataframes
      data.export.ls <- readRDS(paste0(export.data, "TMP.RDS"))
      data.export.ls <- data.export.ls[which(lengths(data.export.ls) != 0)]

      # Model results
      model.export.ls <- readRDS(paste0(export.data, "MOD.RDS"))
      model.export.ls <- model.export.ls[which(lengths(model.export.ls) != 0)]
    }

    # This creates individual dataframes for each experiment
    tmp.ls <- lapply(data.export.ls, function(x){

      result <- x %>%
        group_split(`Experiment Ramp ID`)

      return(result)
    })

    # Flatten lists and remove failed model runs
    {
      data.ls <- flatten(tmp.ls)
      model.ls <- flatten(model.export.ls)

      # Filter failed models
      model.ls <- model.ls[-which(map(model.ls, 6) == "FAIL")]
    }

    # Retrieve dates and experiment numbers
    names.ch <- paste0(as.Date(unlist(map(model.ls, 4))), " ", unlist(map(model.ls, 5)))

    # Rename lists
    names(data.ls) <- names.ch
    names(model.ls) <- names.ch

    # Start loop
    for (n in 1:length(data.ls)){

      # Retrieve data from the list
      dataSPIN.df <- data.ls[[n]]

      # Retrieve model data from the list
      model <- model.ls[[n]]

      # Date string for plotting later
      date.c <- unique(lubridate::as_date(dataSPIN.df$`UTC Time`))[1]

      # Experiment ID
      ID <- unique(dataSPIN.df$`Experiment Ramp ID`)

      bins.ix <- suppressWarnings(which(!is.na(as.numeric(colnames(dataSPIN.df)))))
      bins.nm <- colnames(dataSPIN.df)[bins.ix]

      # Plot path
      export.plot.path = paste0(export.plot, str_remove_all(date.c, '-'), "/")

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Size Histograms #####
      #'

      {
        print(paste0(date.c, ": Size Histograms"))

        title.main <- paste0("Spectrometer for Ice Nucleation (SPIN)")
        title.sub <- paste0("Size Histograms: ", date.c, ", Experiment ", ID)

        # Colorblind accessible color palette
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

        plot.filename <- paste0(export.plot.path, date.c, " Size Histograms", " Experiment ", ID, ".png")

        # Convert to long format for easier plotting
        tmp.long <- dataSPIN.df %>%
          select(`ML Class`, all_of(bins.ix)) %>%
          pivot_longer(
            cols = !c("ML Class"),
            names_to = "Variable",
            values_to = "Value"
          )

        tmp.long <- tmp.long %>%
          mutate(`Variable` = factor(`Variable`, levels = as.character(bins.nm)))

        plot.gg <- ggplot(tmp.long, aes(fill = `ML Class`, y = `Value`, x = `Variable`)) +
          geom_bar(stat = "identity")

        ggsave(
          plot.filename,
          plot.gg,
          width = 8,
          height = 8,
          dpi = 300,
          units = "in",
          bg = "#ffffff"
        )
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Classification Plot #####
      #

      {
        print(paste0(date.c, ": Classification"))

        dataSPIN.df <- dataSPIN.df %>%
          mutate(`Class` = factor(`Class`, levels = c("Aerosol", "Droplet", "Ice", "Water Uptake", "Other")))

        title.main <- paste0("Spectrometer for Ice Nucleation (SPIN)")
        title.sub <- paste0("Classified: ", date.c, ", Experiment ", ID)

        # Colorblind accessible color palette
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

        gg1 <- ggplot(dataSPIN.df, aes(x = `OPC Size`, y = `Depolarization`, col = `Class`)) +
          geom_point(size = 0.1, alpha = 0.5) +
          scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          scale_color_manual(values = cbPalette) +
          labs(title = "Preliminary Classification") +
          theme(
            plot.title = element_text(size = 10, hjust = 0.5),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_text(vjust = -1.5, size = 10),
            axis.title.y = element_text(vjust = 2, size = 10),
            axis.ticks.x = element_line(linewidth = 0.5),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 0.2, 0.25, 0.2), "cm"),
            aspect.ratio = 1
          ) + coord_cartesian(clip = "off") +
          guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

        gg2 <- ggplot(dataSPIN.df, aes(x = `OPC Size`, y = `Depolarization`, col = `ML Class`)) +
          geom_point(size = 0.1, alpha = 0.5) +
          scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          scale_color_manual(values = cbPalette) +
          labs(title = "PCA-SVM Classification") +
          theme(
            plot.title = element_text(size = 10, hjust = 0.5),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_text(vjust = -1.5, size = 10),
            axis.title.y = element_text(vjust = 2, size = 10),
            axis.ticks.x = element_line(linewidth = 0.5),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 0.2, 0.25, 0.2), "cm"),
            aspect.ratio = 1
          ) + coord_cartesian(clip = "off") +
          guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

        gg3 <- ggplot(dataSPIN.df, aes(x = `Lamina S Ice`, y = `Depolarization`, col = `Class`)) +
          geom_point(size = 0.1, alpha = 0.5) +
          scale_x_continuous(limits = c(1, 1.7), breaks = seq(1, 1.7, 0.1)) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          scale_color_manual(values = cbPalette) +
          theme(
            plot.title = element_text(),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_text(vjust = -1.5, size = 10),
            axis.title.y = element_text(vjust = 2, size = 10),
            axis.ticks.x = element_line(linewidth = 0.5),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 0.2, 0.25, 0.2), "cm"),
            aspect.ratio = 1
          ) + coord_cartesian(clip = "off") +
          guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

        gg4 <- ggplot(dataSPIN.df, aes(x = `Lamina S Ice`, y = `Depolarization`, col = `ML Class`)) +
          geom_point(size = 0.1, alpha = 0.5) +
          scale_x_continuous(limits = c(1, 1.7), breaks = seq(1, 1.7, 0.1)) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          scale_color_manual(values = cbPalette) +
          theme(
            plot.title = element_text(),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_text(vjust = -1.5, size = 10),
            axis.title.y = element_text(vjust = 2, size = 10),
            axis.ticks.x = element_line(linewidth = 0.5),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 0.2, 0.25, 0.2), "cm"),
            aspect.ratio = 1
          ) + coord_cartesian(clip = "off") +
          guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

        total.gg <- gg1 + gg2 + gg3 + gg4 + plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") +
          plot_annotation(
            title = paste0('Classification Performance - ', date.c)
          )

        plot.filename <- paste0(export.plot, date.c, " Classification.png")

        ggsave(filename = plot.filename, total.gg, width = 10, height = 8)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Model Performance #####
      #'

      {
        # Model Performance
        {
          print(paste0(date.c, ": ML Model Performance"))

          plot.filename <- paste0(export.plot.path, date.c, " ML Model Performance", " Experiment", ".png")

          tmp <- ggplot(model[[1]]$results, aes(x = C, y = Accuracy)) +
            geom_point() +
            geom_line(color = 'blue') +
            theme_minimal()

          ggsave(filename = plot.filename, tmp, width = 8, height = 6, bg = "white")
        }

        # Model Performance
        {
          print(paste0(date.c, ": PCA Model Performance"))

          plot.filename <- paste0(export.plot.path, date.c, " PCA Model Performance", " Experiment", ".png")

          tmp <- factoextra::fviz_pca_biplot(model[[3]]) + factoextra::fviz_screeplot(model[[3]])

          ggsave(filename = plot.filename, tmp, width = 8, height = 6)
        }
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Light Scattering  #####
      #'

      {
        tmp.df <- dataSPIN.df %>%
          filter(`Inlet Filter ON` == 0)

        plot.filename <- paste0(export.plot.path, date.c, " Scattering Classification", " Experiment ", ID, ".png")

        annotation.df <- data.frame(Temperature = unique(as.numeric(as.character(tmp.df$`Lamina Breaks`))))

        annotation.df <- annotation.df %>%
          mutate(`Annotation` = homogeneous.freezing.solver(Temperature, 0.3, 60)) %>%
          mutate(`Annotation` = if_else(`Temperature` > -38, (100*p_liquid(Temperature)/p_ice(Temperature))/100, Annotation)) %>%
          mutate(`Label` = if_else(`Temperature` > -38, "Water Saturation", "Homogeneous Freezing"))

        gg1 <- ggplot(tmp.df, aes(y = `Depolarization`, x = after_stat(scaled), fill = as.factor(`Lamina Breaks`), color = as.factor(`Lamina Breaks`), group = `Lamina Breaks`)) +
          geom_density(alpha = 0.125, adjust = 2.5) +
          scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
          scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
          ylab(label = latex2exp::TeX(r'($\delta_{SPIN}$)')) +
          xlab(label = "Kernel Density Estimate (KDE)") +
          scale_fill_manual(values = c("#E69F00", "#999999", "#56B4E9", "#009E73")) +
          scale_color_manual(values = c("#E69F00", "#999999", "#56B4E9", "#009E73")) +
          annotate("text", y = 0.9, x = 0.9, label = "A", fontface = "bold", size = 6) +
          theme(panel.background = element_rect(fill = "white"),
                plot.background = element_rect(fill = "white"),
                panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
                panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
                panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
                axis.title.x = element_text(vjust = -1, size = 12),
                axis.title.y = element_text(vjust = 2, size = 12),
                axis.ticks.length = unit(2, "mm"),
                axis.minor.ticks.length = unit(1, "mm"),
                legend.title = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
          coord_fixed()

        gg2 <- ggplot(tmp.df, aes(x = `Lamina S Ice`, y = `Depolarization`, color = as.factor(`Lamina Breaks`), group = `Lamina Breaks`)) +
          annotate("rect", ymin = 0.4, ymax = 1, xmin = 1, xmax = 1.6, fill = "gray90", alpha = 0.2) +
          annotate("rect", ymin = 0.16, ymax = 0.4, xmin = 1, xmax = 1.6, fill = "gray60", alpha = 0.2) +
          annotate("rect", ymin = 0, ymax = 0.2, xmin = 1, xmax = 1.6, fill = "white", alpha = 0.2) +
          geom_point(size = 0.01, alpha = 0.1) +
          geom_smooth(alpha = 0.5) +
          geom_vline(data = annotation.df, aes(xintercept = `Annotation`, col = as.factor(`Temperature`)), linetype = 2) +
          geom_text(data = annotation.df, aes(x = `Annotation`, y = 0.8, label = `Label`, color = as.factor(`Temperature`)), inherit.aes = F,
                    angle = 90, nudge_x = -0.0105, size = 3) +
          scale_x_continuous(breaks = seq(1, 1.6, 0.1), limits = c(1, 1.6), expand = c(0.001, 0.001)) +
          scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
          xlab(label = latex2exp::TeX(r'($S_{ice}$)')) +
          scale_color_manual(values = c("#E69F00", "#999999", "#56B4E9", "#009E73")) +
          annotate("text", y = 0.9, x = 1.55, label = "B", fontface = "bold", size = 6) +
          annotate("text", y = 0.55, x = 1.025, label = "Ice", hjust = 0, size = 8/.pt) +
          annotate("text", y = 0.35, x = 1.025, label = "Droplet", hjust = 0, size = 8/.pt) +
          annotate("text", y = 0.11, x = 1.025, label = "Aerosol/Water Uptake", hjust = 0, size = 8/.pt) +
          theme(panel.background = element_rect(fill = "white"),
                plot.background = element_rect(fill = "white"),
                panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
                panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
                panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
                axis.text.y = element_blank(),
                axis.title.x = element_text(vjust = -2, size = 12),
                axis.title.y = element_blank(),
                axis.ticks.length = unit(2, "mm"),
                axis.minor.ticks.length = unit(1, "mm"),
                panel.ontop = F,
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
          guides(color = "none") +
          coord_fixed(ratio = 1/1.6)

        plot.gg <- gg1 + gg2 + plot_layout(guides = "collect") & theme(legend.position = "right")

        ggsave(
          plot.filename,
          plot.gg,
          width = 10,
          height = 5,
          dpi = 300,
          units = "in",
          bg = "#ffffff"
        )
      }
    }
  }
}
