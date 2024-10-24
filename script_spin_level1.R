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

  # Plotting libraries
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(ggsci)
  library(viridis)
  library(RColorBrewer)

  # Plotting parameters
  resolution.dpi = 400
  font.family = "Helvetica"

  # Set working directory
  setwd("~/Library/CloudStorage/Box-Box/Organosulfate Proxies")
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
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_classifier.R"))
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

{
  # Plotting on?
  PLOT.ON = T

  # Minimum threshold in um to be considered ice
  THRESHOLD.Dp.nm = 2.5
  {
    data.export.ls <- NULL
    for (n in 1:length(spin.dates)){

      {
        print(paste0("Level 1 SPIN Data for ", spin.dates[n], " from directory ", spin.path[n]))

        # Memory management
        gc(verbose = F)

        # Code Timing
        pct <- proc.time()

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
                mutate(`Total Particles` = rowSums(select(., .dots = all_of(cutoff.bins.nm)))) %>%
                mutate(`Inlet Filter ON` = dataSPIN.df$`Inlet Filter ON`) %>%
                select(`Inlet Filter ON`, `Total Particles`)

              # Set total backgrounds to NA for when the inlet filter is on
              # Using 0 will sway statistics significantly
              tmp.df <- tmp.df %>%
                mutate(`Total Background` = if_else(`Inlet Filter ON` == 1, `Total Particles`, NA), .after = `Total Particles`)

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
                select(!`Inlet Filter ON`)

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
              relocate(c(`Total Particles`, `Total Background`), .after = `Units`)
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

            # Find lamina breaks that segment dataset
            dataALL.df <- dataALL.df %>%
              mutate(`Lamina Breaks` = round(`Lamina Temp (C)`/50, 1)*50)

            # Calculate the depolarization ratio
            # Normalize values to be between 0 and 1 for standardization
            dataALL.df <- dataALL.df %>%
              mutate(`Depolarization` = `S1`/(`P1` + `P2`)) %>%
              mutate(`Depolarization` = if_else(Depolarization > 1, NA, Depolarization)) %>%
              mutate(`Depolarization` = if_else(Depolarization < 0, NA, Depolarization)) %>%
              filter(!is.na(`Depolarization`))

            # Order columns
            dataALL.df <- dataALL.df %>%
              relocate(`S1`:`Depolarization`, .after = `Inlet Filter ON`)

            # Date string for plotting later
            date.c <- unique(lubridate::as_date(dataALL.df$`UTC Time`))[1]

            print(paste0(date.c, ": Data Aggregation Complete"))

            # Clean up memory
            rm(dataBINS.df, rawPBP.df, rawSPIN.df)
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

          {
            print(paste0(date.c, ": Classification"))

            # Apply classifier to dataset
            class.df <- spin.classifier(dataALL.df, ML.option = F, size.limit = THRESHOLD.Dp.nm)

            # Reorder columns
            tmp.c <- colnames(class.df)

            # apply classifier function
            export.df <- cbind(dataALL.df, class.df) %>%
              relocate(all_of(tmp.c), .after = "Depolarization")
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
      }
    }

    # Backup data from current run for testing
    saveRDS(data.export.ls, file = paste0(export.data, "TMP.RDS"))

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

    # Read back in looped data from above and filter empty data
    data.export.ls <- readRDS(paste0(export.data, "TMP.RDS"))
    data.export.ls <- data.export.ls[which(lengths(data.export.ls) != 0)]

    # This creates individual dataframes for each experiment
    data.export.ls <- lapply(data.export.ls, function(x){

      result <- x %>%
        group_split(`Experiment Ramp ID`)

      return(result)
    })

    # Flatten the lists to remove nested levels
    data.export.ls <- purrr::list_flatten(data.export.ls)

    # Start loop
    for (n in 1:length(data.export.ls)){

      # Retrieve data from the list
      dataSPIN.df <- data.export.ls[[n]]

      # Date string for plotting later
      date.c <- unique(lubridate::as_date(dataSPIN.df$`UTC Time`))[1]

      # Experiment ID
      ID <- unique(dataSPIN.df$`Experiment Ramp ID`)

      bins.ix <- suppressWarnings(which(!is.na(as.numeric(colnames(dataSPIN.df)))))
      bins.nm <- colnames(dataSPIN.df)[bins.ix]

      # Plot path
      export.plot.path = paste0(export.plot, str_remove_all(date.c, '-'), "/")

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: Classification Plot #####
      #

      {
        print(paste0(date.c, ": Classification"))

        title.main <- paste0("Spectrometer for Ice Nucleation (SPIN)")
        title.sub <- paste0("Classified: ", date.c, ", Experiment ", ID)

        # Colorblind accessible color palette
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

        plot.filename <- paste0(export.plot.path, date.c, " Classification", " Experiment ", ID, ".png")

        gg1 <- ggplot(dataSPIN.df, aes(x = `Lamina S Ice`, y = `Depolarization`, color = `Class`)) +
          geom_point(size = 0.1, alpha = 0.5) +
          scale_x_continuous(limits = c(1, 1.7)) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          xlab(latex2exp::TeX(r'(Lamina $S_{ice}$)')) +
          scale_color_manual(values = cbPalette) +
          facet_wrap(~`Lamina Breaks`) +
          theme(
            plot.title = element_text(),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_text(vjust = -1.5),
            axis.title.y = element_text(vjust = 2),
            axis.ticks.x = element_line(linewidth = 0.5),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 0.2, 0.25, 0.2), "cm"),
            panel.spacing = unit(2, "lines")
          ) + coord_cartesian(clip = "off") +
          guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

        gg2 <- ggplot(dataSPIN.df, aes(x = `Log Size`, y = `Depolarization`, color = `Class`)) +
          geom_point(size = 0.1, alpha = 0.5) +
          scale_x_continuous(limits = c(1, 5), expand = c(0, 0)) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
          facet_wrap(~`Lamina Breaks`) +
          scale_color_manual(values = cbPalette) +
          theme(
            plot.title = element_text(),
            plot.subtitle = element_text(color = "gray25"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
            panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
            panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_text(vjust = -1.5),
            axis.title.y = element_text(vjust = 2),
            axis.ticks.x = element_line(linewidth = 0.5),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.length = unit(1, "mm"),
            plot.margin = unit(c(0.2, 0.2, 0.25, 0.2), "cm"),
            panel.spacing = unit(2, "lines")
          ) + coord_cartesian(clip = "off") +
          guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

        plot.gg <- ggarrange(gg1, gg2, nrow = 2, common.legend = T, legend = "right") +
          plot_annotation(title = title.main, subtitle = title.sub)

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

        gg1 <- ggplot(tmp.df, aes(y = `Depolarization`, x = after_stat(scaled), fill = as.factor(`Lamina Breaks`), group = `Lamina Breaks`)) +
          geom_density(alpha = 0.4, adjust = 2.5) +
          scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
          scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
          ylab(label = latex2exp::TeX(r'($\delta_{SPIN}$)')) +
          xlab(label = "Kernel Density Estimate (KDE)") +
          scale_fill_manual(values = cbPalette) +
          theme(panel.background = element_rect(fill = "white"),
                plot.background = element_rect(fill = "white"),
                panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
                panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
                panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
                axis.title.x = element_text(vjust = -1),
                axis.title.y = element_text(vjust = 2, size = 14),
                axis.ticks.length = unit(2, "mm"),
                axis.minor.ticks.length = unit(1, "mm"),
                legend.title = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
          coord_fixed()

        gg2 <- ggplot(tmp.df, aes(x = `Lamina S Ice`, y = `Depolarization`, color = as.factor(`Lamina Breaks`), group = `Lamina Breaks`)) +
          annotate("rect", ymin = 0.4, ymax = 0.6, xmin = 1, xmax = 1.6, fill = "gray60", alpha = 0.2) +
          annotate("text", y = 0.5, x = 1.05, label = "Ice") +
          annotate("rect", ymin = 0.15, ymax = 0.4, xmin = 1, xmax = 1.6, fill = "gray80", alpha = 0.2) +
          annotate("text", y = 0.275, x = 1.05, label = "Droplet") +
          geom_point(size = 0.01, alpha = 0.1) +
          geom_smooth(alpha = 0.5) +
          geom_vline(data = annotation.df, aes(xintercept = `Annotation`, col = as.factor(`Temperature`)), linetype = 2) +
          geom_text(data = annotation.df, aes(x = `Annotation`, y = 0.8, label = `Label`, color = as.factor(`Temperature`)), inherit.aes = F,
                    angle = 90, nudge_x = -0.0105, size = 3) +
          scale_x_continuous(breaks = seq(1, 1.6, 0.1), limits = c(1, 1.6), expand = c(0.001, 0.001)) +
          scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
          xlab(label = latex2exp::TeX(r'($S_{ice}$)')) +
          scale_color_manual(values = cbPalette) +
          theme(panel.background = element_rect(fill = "white"),
                plot.background = element_rect(fill = "white"),
                panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
                panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
                panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
                axis.text.y = element_blank(),
                axis.title.x = element_text(vjust = -2),
                axis.title.y = element_blank(),
                axis.ticks.length = unit(2, "mm"),
                axis.minor.ticks.length = unit(1, "mm"),
                panel.ontop = F,
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
          guides(color = "none") +
          coord_fixed(ratio = 1/1.6)

        plot.gg <- (gg1 & theme(plot.tag.position  = c(.935, .96))) -
          (gg2 & theme(plot.tag.position  = c(.785, .96))) +
          plot_annotation(tag_levels = "a") +
          plot_layout(guides = "collect") & theme(legend.position = "right")

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

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: Correlation Plots #####
      #'

      {
        print(paste0(date.c, ": Plotting Correlation Plots"))

        variables.c <- c(bins.nm, "Lamina S Liquid", "Lamina S Ice",
                         "Log S1", "Log P1", "Log P2", "Log Size", "Depolarization")

        X1.df <- dataSPIN.df %>%
          filter(`Class` == "Aerosol") %>%
          select(all_of(variables.c))

        X2.df <- dataSPIN.df %>%
          filter(`Class` == "Ice") %>%
          select(all_of(variables.c))

        X3.df <- dataSPIN.df %>%
          filter(`Class` == "Droplet") %>%
          select(all_of(variables.c))

        X4.df <- dataSPIN.df %>%
          filter(`Class` == "Other") %>%
          select(all_of(variables.c))

        C1 <- cor(X1.df, use = "complete.obs")
        C2 <- cor(X2.df, use = "complete.obs")
        C3 <- cor(X3.df, use = "complete.obs")
        C4 <- cor(X4.df, use = "complete.obs")

        plot.filename = paste0(export.plot.path, date.c, " Correlation Matrices", " Experiment ", ID, ".pdf")

        pdf(
          file = plot.filename,
          width = 10,
          height = 10,
          bg = "white"
        )

        margins = c(0.1, 1, 1, 1)

        {
          title.main = paste0("SPIN Correlation Matrix: ", date.c, " Experiment ", ID)
          title.sub = "Aerosol"

          corrplot::corrplot(C1, method = 'circle', cl.align.text = "l", type = "full", order = "original", diag = TRUE,
                             tl.col = "black", tl.srt = 90, tl.cex = 1, cl.cex = 1, cl.offset = 0.1,
                             tl.offset = 0.3, insig = "pch", mar = margins,
                             na.label = "o", na.label.col = "black")

          mtext(side = 3, line = 0, at = -0.07, adj = 0, cex = 1.5, title.main, family = 'Helvetica')
          mtext(side = 3, line = -1.5, at = -0.07, adj = 0, cex = 1, title.sub)
        }

        {
          title.main = paste0("SPIN Correlation Matrix: ", date.c, " Experiment ", ID)
          title.sub = "Ice"

          corrplot::corrplot(C2, method = 'circle', cl.align.text = "l", type = "full", order = "original", diag = TRUE,
                             tl.col = "black", tl.srt = 90, tl.cex = 1, cl.cex = 1, cl.offset = 0.1,
                             tl.offset = 0.3, insig = "pch", mar = margins,
                             na.label = "o", na.label.col = "black")

          mtext(side = 3, line = 0, at = -0.07, adj = 0, cex = 1.5, title.main, family = 'Helvetica')
          mtext(side = 3, line = -1.5, at = -0.07, adj = 0, cex = 1, title.sub)
        }

        {
          title.main = paste0("SPIN Correlation Matrix: ", date.c, " Experiment ", ID)
          title.sub = "Droplet"

          corrplot::corrplot(C3, method = 'circle', cl.align.text = "l", type = "full", order = "original", diag = TRUE,
                             tl.col = "black", tl.srt = 90, tl.cex = 1, cl.cex = 1, cl.offset = 0.1,
                             tl.offset = 0.3, insig = "pch", mar = margins,
                             na.label = "o", na.label.col = "black")

          mtext(side = 3, line = 0, at = -0.07, adj = 0, cex = 1.5, title.main, family = 'Helvetica')
          mtext(side = 3, line = -1.5, at = -0.07, adj = 0, cex = 1, title.sub)
        }

        {
          title.main = paste0("SPIN Correlation Matrix: ", date.c, " Experiment ", ID)
          title.sub = "Other"

          corrplot::corrplot(C4, method = 'circle', cl.align.text = "l", type = "full", order = "original", diag = TRUE,
                             tl.col = "black", tl.srt = 90, tl.cex = 1, cl.cex = 1, cl.offset = 0.1,
                             tl.offset = 0.3, insig = "pch", mar = margins,
                             na.label = "o", na.label.col = "black")

          mtext(side = 3, line = 0, at = -0.07, adj = 0, cex = 1.5, title.main, family = 'Helvetica')
          mtext(side = 3, line = -1.5, at = -0.07, adj = 0, cex = 1, title.sub)
        }

        dev.off()
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: Probability Density of Uni-variate PBP Data #####
      #'

      {
        print(paste0(date.c, ": Plotting Probability Density of Uni-variate PBP Data"))

        {
          D1.df <- X1.df %>%
            mutate(`Index` = seq(1, n(), 1)) %>%
            select("Index", "Log S1", "Log P1", "Log P2", "Log Size")

          D2.df <- X2.df %>%
            mutate(`Index` = seq(1, n(), 1)) %>%
            select("Index", "Log S1", "Log P1", "Log P2", "Log Size")

          D3.df <- X3.df %>%
            mutate(`Index` = seq(1, n(), 1)) %>%
            select("Index", "Log S1", "Log P1", "Log P2", "Log Size")

          D4.df <- X4.df %>%
            mutate(`Index` = seq(1, n(), 1)) %>%
            select("Index", "Log S1", "Log P1", "Log P2", "Log Size")
        }

        x.breaks = pretty(c(-1, 5), 5)
        x.lims = range(x.breaks)

        plot.filename = paste0(export.plot.path, date.c, " PDFs", " Experiment ", ID, ".pdf")

        pdf(
          file = plot.filename,
          width = 10,
          height = 10,
          bg = "white"
        )

        plot.list <- list()
        {
          title.main = paste0("Optical Particle Counter - Probability Density Functions: ", date.c, " Experiment ", ID)
          title.sub = "Aerosol"

          plot.list[[1]] <- density.histogram.png(title.main,
                                                  title.sub,
                                                  data = D1.df,
                                                  index.vars = "Index",
                                                  x.breaks,
                                                  x.lims,
                                                  plot.width = 10,
                                                  plot.height = 8,
                                                  plot.resolution.dpi = 400)
        }

        {
          title.main = paste0("Optical Particle Counter - Probability Density Functions: ", date.c, " Experiment ", ID)
          title.sub = "Ice"

          plot.list[[2]] <- density.histogram.png(title.main,
                                                  title.sub,
                                                  data = D2.df,
                                                  index.vars = "Index",
                                                  x.breaks,
                                                  x.lims,
                                                  plot.width = 10,
                                                  plot.height = 8,
                                                  plot.resolution.dpi = 400)
        }

        {
          title.main = paste0("Optical Particle Counter - Probability Density Functions: ", date.c, " Experiment ", ID)
          title.sub = "Droplet"

          plot.list[[3]] <- density.histogram.png(title.main,
                                                  title.sub,
                                                  data = D3.df,
                                                  index.vars = "Index",
                                                  x.breaks,
                                                  x.lims,
                                                  plot.width = 10,
                                                  plot.height = 8,
                                                  plot.resolution.dpi = 400)
        }

        {
          title.main = paste0("Optical Particle Counter - Probability Density Functions: ", date.c, " Experiment ", ID)
          title.sub = "Other"

          plot.list[[4]] <- density.histogram.png(title.main,
                                                  title.sub,
                                                  data = D4.df,
                                                  index.vars = "Index",
                                                  x.breaks,
                                                  x.lims,
                                                  plot.width = 10,
                                                  plot.height = 8,
                                                  plot.resolution.dpi = 400)
        }

        invisible(lapply(plot.list, print))
        dev.off()
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: Bivariate KDE Joint Probability Density #####
      #'
      #'

      {
        print(paste0(date.c, ": Plotting Bivariate KDE Joint Probability Density"))

        {
          title.main = ""
          title.sub = "Aerosol"

          gg1 <- density.contour.png(title.main,
                                     title.sub,
                                     data = X1.df,
                                     binwidth.c = 0.05,
                                     raster.n = 200,
                                     plot.filename = NULL,
                                     plot.width = 8,
                                     plot.height = 8,
                                     plot.resolution.dpi = 400)
        }

        {
          title.main = ""
          title.sub = "Ice"

          gg2 <- density.contour.png(title.main,
                                     title.sub,
                                     data = X2.df,
                                     binwidth.c = 0.05,
                                     raster.n = 200,
                                     plot.filename = NULL,
                                     plot.width = 8,
                                     plot.height = 8,
                                     plot.resolution.dpi = 400)
        }

        {
          title.main = ""
          title.sub = "Droplet"

          gg3 <- density.contour.png(title.main,
                                     title.sub,
                                     data = X3.df,
                                     binwidth.c = 0.05,
                                     raster.n = 200,
                                     plot.filename = NULL,
                                     plot.width = 8,
                                     plot.height = 8,
                                     plot.resolution.dpi = 400)
        }

        {
          title.main = ""
          title.sub = "Other"

          gg4 <- density.contour.png(title.main,
                                     title.sub,
                                     data = X4.df,
                                     binwidth.c = 0.05,
                                     raster.n = 200,
                                     plot.filename = NULL,
                                     plot.width = 8,
                                     plot.height = 8,
                                     plot.resolution.dpi = 400)
        }

        plot.filename = paste0(export.plot.path, date.c, " OPC Joint PDFs", " Experiment ", ID, ".png")

        plot.gg <- ggpubr::ggarrange(gg1 + rremove("ylab") + rremove("xlab"),
                                     gg4 + rremove("ylab") + rremove("xlab"),
                                     gg3 + rremove("ylab") + rremove("xlab"),
                                     gg2 + rremove("ylab") + rremove("xlab"),
                                     nrow = 2, ncol = 2,
                                     widths = c(1, 1), common.legend = T, legend = "right",
                                     labels = NULL)

        plot.gg <-
          annotate_figure(plot.gg,
                          top = text_grob(paste0("Joint Probability Densities - 2D Gaussian KDE: ", date.c, " Experiment ", ID), size = 14),
                          left = text_grob(latex2exp::TeX(r'($Depolarization$)', bold = F), rot = 90, size = 10),
                          bottom = text_grob(latex2exp::TeX(r'($Log_{10}(Size)$)', bold = F), rot = 0, size = 10))

        ggsave(
          plot.filename,
          plot.gg,
          width = 12,
          height = 10,
          dpi = 400,
          units = "in",
          bg = "#ffffff"
        )
      }
    }
  }
}
