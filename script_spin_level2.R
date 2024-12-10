##### SCRIPT: SPectrometer for Ice Nucleation (SPIN) Level 2 #####
#' @author Christopher Rapp
#' @description
#'
#' NOTE: This code has been optimized with parallel processing. If you are
#' struggling with errors please swap to a normal for-loop or lapply rather than mcapply
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
  library(patchwork)

  # These libraries are needing to implement parallel programming in R
  library(parallel)
  library(foreach)
  library(doParallel)

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
  import.spin = paste0(work.dir, "/SPIN/export/level1/")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/SPIN/export/level2/")
  export.plot = paste0(work.dir, "/SPIN/export/plots/")

  #' @importFrom
  #' Function scripts
  source(paste0("~/Documents/GitHub/functions/", "functions_aerosols.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_microphysics.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_plotting.R"))
  source(paste0("~/Documents/GitHub/functions/", "functions_spin_lamina.R"))
}

# ---------------------------------------------------------------------------- #
##### SECTION: Parallel Computing Setup #####
#'

{
  # Use 75% of the available cores of the machine running this code
  CORES.available <- detectCores()
  CORES.dedicated = as.integer(floor(CORES.available*0.75))
}


# ---------------------------------------------------------------------------- #
##### SECTION: File Selection #####
#'

{
  # Directories
  spin.files <- list.files(path = import.spin, recursive = F, full.names = T, pattern = ".csv")

  # Pull POSIX style dates from directory names
  spin.dates <- as_date(str_extract(spin.files, '(?<=SPIN_)\\d{8}'), format = '%Y%m%d')
}

# ---------------------------------------------------------------------------- #
##### SECTION: Level 2 #####
#'

{
  # ---------------------------------------------------------------------- #
  ##### SPIN Data #####
  #

  {
    print(paste0("SPIN Data"))

    dataSPIN.ls <- lapply(spin.files, function(x){

      # Loop through SPIN files and read in
      # SPIN is exported using f.write which preserves time zone in the file
      rawSPIN.df <- fread(x)

      # Filter data after an experimental ramp has started
      rawSPIN.df <- rawSPIN.df %>%
        mutate(`Local Time` = lubridate::with_tz(`Local Time`, tzone = tz.c)) %>%
        mutate(`Start Time` = lubridate::with_tz(`Local Time`, tzone = tz.c))

      # Set lamina breaks to a factor
      dataSPIN.df <- rawSPIN.df %>%
        mutate(`Lamina Breaks` = as.factor(`Lamina Breaks`)) %>%
        relocate(c(`Lamina Breaks`, `Depolarization`), .after = `Lamina S Liquid`) %>%
        relocate(`Lamina Temp (C)`, .after = `Lamina S Liquid`)

      return(dataSPIN.df)
    })

    # Split data into experimental groups
    dataSPIN.ls <- lapply(dataSPIN.ls, function(x){
      split(x, x$`Experiment Ramp ID`, lex.order = T)
    })

    # Flatten list
    dataSPIN.ls <- purrr::list_flatten(dataSPIN.ls)

    # Rename list elements
    dataSPIN.ls <- setNames(dataSPIN.ls, nm = seq(1, length(dataSPIN.ls), 1))
  }

  # -------------------------------------------------------------------------- #
  ##### SECTION: Lamina Uncertainty #####
  #

  {
    print(paste0("Lamina Error Analysis"))

    dataSPIN.ls <- mclapply(dataSPIN.ls, mc.cores = CORES.dedicated, function(x){

      # Calculate error propagation
      error.df <- lamina.error(x, position.error = T)

      # Add calculated error propagation values
      data.df <- cbind(x, error.df) %>%
        relocate(all_of(colnames(error.df)), .after = `Inlet Filter ON`)

      # Remove error points when there is an issue with the warm wall
      data.df <- data.df %>%
        mutate(`Lamina S Ice Error` = if_else(`Inlet Filter ON` == 1, NA, `Lamina S Ice Error`))

      return(data.df)
    })
  }

  # -------------------------------------------------------------------------- #
  ##### SECTION: Chamber Analysis #####
  #' Thermocouple pair by pair comparison of expected lamina conditions
  #

  print(paste0("Chamber Analysis"))

  dataSPIN.ls <- mclapply(dataSPIN.ls, mc.cores = CORES.dedicated, function(x){

    # Select spin.dates from loop
    date.c <- unique(as.Date(x$`Local Time`))

    # Find unique experiment ramp
    ID.c <- unique(x$`Experiment Ramp ID`)

    # Export plot path
    export.plot.path = paste0(export.plot, stringr::str_remove_all(date.c, "-"), '/')

    # Smallest diameter that could be a heterogeneous ice nuclei from Vali
    r_um = 0.1

    # Run lamina conditions function on the chamber wall thermocouples.
    # Outputs a list of length 3
    # List 1 - Lamina S ice
    # List 2 - Lamina S liquid
    # List 3 - Homogeneous Freezing Thresholds
    tmp.ls <- lamina.conditions(x, x$`Lamina Centerline Ratio`, r_um = 0.1, t = 60)

    # Find difference in predicted lamina s ice compared to each thermocouple
    tmp.df <- cbind("UTC Time" = tmp.ls[[1]][, 1], tmp.ls[[1]][, 2:17] - x$`Lamina S Ice`)

    caption1.c = paste0("SPIN chamber conditions using paired cold-warm wall thermocouples
                          i.e. C7 & W7. Conditions are generated using a linear interpolation of
                          saturation vapor pressures and temperatures per thermocouple for each
                          experiment time step. Vapor pressure calculations performed as described by Murphy and Koop (2005).
                          Distances assumed to be 10 mm and neglecting change due to ice surface.
                          Values are selected from interpolation vector by using the
                          calculated lamina position as extended by Rogers (1988). This position is
                          designated specifically for centerline thermocouples
                          T0, T3, T6, T9, T12 but extended to the entire chamber to locate any
                          inhomogeneities. Gray hollow circles mark the expected homogeneous freezing
                          threshold as detailed by Koop (2002) for ", r_um ," \U00B5m particles at the
                          linearly interpolated lamina temperature.")

    caption2.c = paste0("SPIN chamber conditions using paired cold-warm wall thermocouples
                          i.e. C7 & W7. Conditions are generated using a linear interpolation of
                          saturation vapor pressures and temperatures per thermocouple for each
                          experiment time step. Vapor pressure calculations performed as described by Murphy and Koop (2005).
                          Distances assumed to be 10 mm and neglecting change due to ice surface.
                          Values are selected from interpolation vector by using the
                          calculated lamina position as extended by Rogers (1988). This position is
                          designated specifically for centerline thermocouples
                          T0, T3, T6, T9, T12 but extended to the entire chamber to locate any
                          inhomogeneities.")

    caption3.c = paste0("Difference between ice supersaturations derived from
                           individually paired cold-warm wall thermocouples
                           and reported SPIN lamina conditions. Lower supersaturation
                           values at the bottom of the chamber e.g. T12 - T15 are due
                           to the droplet evaporation section.")

    # Chamber ice saturation
    chamber.contour.png(data = tmp.ls[[1]],
                        data2 = tmp.ls[[3]],
                        title = paste0("Spectrometer for Ice Nucleation (SPIN)"),
                        subtitle = paste0("Chamber Ice Saturation: ", date.c),
                        title.key = latex2exp::TeX(r'($S_{ice}$)', bold = T),
                        levels.limits = c(1, 2),
                        col.palette = NULL,
                        n.colors = 20,
                        plot.filename = paste0(export.plot.path, date.c, ' Chamber Ice Saturation Experiment ', ID.c, '.png'),
                        plot.width = 12,
                        plot.height = 6,
                        plot.resolution.dpi = 600,
                        caption = caption1.c)

    # Chamber liquid saturation
    chamber.contour.png(data = tmp.ls[[2]],
                        data2 = NULL,
                        title = paste0("Spectrometer for Ice Nucleation (SPIN)"),
                        subtitle = paste0("Chamber Liquid Saturation: ", date.c),
                        title.key = latex2exp::TeX(r'($S_{liq}$)', bold = T),
                        levels.limits = c(0.5, 1.5),
                        n.colors = 20,
                        col.palette = NULL,
                        plot.filename = paste0(export.plot.path, date.c, ' Chamber Liquid Saturation Experiment ', ID.c, '.png'),
                        plot.width = 12,
                        plot.height = 6,
                        plot.resolution.dpi = 600,
                        caption = caption2.c)

    # Chamber Observed Ice Saturation Difference
    chamber.contour.png(data = tmp.df,
                        data2 = NULL,
                        title = paste0("Spectrometer for Ice Nucleation (SPIN)"),
                        subtitle = paste0("Chamber Observed Ice Saturation Difference: ", date.c),
                        title.key = latex2exp::TeX(r'($\Delta S_{ice}$)', bold = T),
                        levels.limits = c(-0.4, 0.4),
                        col.palette = c(rev(RColorBrewer::brewer.pal(7, 'Reds')), "white", "white", RColorBrewer::brewer.pal(7, 'Blues')),
                        n.colors = NULL,
                        plot.filename = paste0(export.plot.path, date.c, ' Chamber Ice Difference Experiment ', ID.c, '.png'),
                        plot.width = 12,
                        plot.height = 6,
                        plot.resolution.dpi = 600,
                        caption = caption3.c)

    # ICE
    {
      # Only keep the first 13 thermcouples
      tmp.df <- tmp.ls[[1]][, c(2:14)]

      # Change names from simple numeric labels
      setnames(tmp.df, new = paste0("S Ice Pair ", colnames(tmp.df)))

      # Find the maximum S for each row
      tmp.df$`Maximum S Ice` <- apply(tmp.df, 1, max)

      # Merge
      data.df <- cbind(x, tmp.df) %>%
        relocate(all_of(colnames(tmp.df)), .before = `Timestamp`) %>%
        relocate(`Maximum S Ice`, .before = `Lamina Position`)
    }

    # LIQUID
    {
      # Only keep the first 13 thermcouples
      tmp.df <- tmp.ls[[2]][, c(2:14)]

      # Change names from simple numeric labels
      setnames(tmp.df, new = paste0("S Liq Pair ", colnames(tmp.df)))

      # Find the maximum S for each row
      tmp.df$`Maximum S Liq` <- apply(tmp.df, 1, max)

      # Merge
      data.df <- cbind(data.df, tmp.df) %>%
        relocate(all_of(colnames(tmp.df)), .before = `Timestamp`) %>%
        relocate(`Maximum S Liq`, .before = `Lamina Position`)
    }

    return(data.df)
  })

  # -------------------------------------------------------------------------- #
  ##### SECTION: Export Level 2 Data #####

  {
    for (i in 1:length(dataSPIN.ls)){

      # Extract loop element
      export.df <- dataSPIN.ls[[i]]

      # Select spin.dates from loop
      date.c <- unique(as.Date(export.df$`Local Time`))

      # Find unique experiment ramp
      ID.c <- unique(export.df$`Experiment Ramp ID`)

      # Check if export path exists
      # If it does not, create it
      if (!dir.exists(export.data)) {
        # Create a dated directory to send plots to
        dir.create(export.data, mode = "777")
      }

      export.filename = paste0(export.data, "SPIN001_SPIN_", stringr::str_remove_all(date.c, '-'), "_Experiment", ID.c, "_level2.csv")

      print(paste0("Exporting Level 2 Data: ", export.filename))

      # Save data using data.table::fwrite
      data.table::fwrite(export.df, file = export.filename, showProgress = T)
    }

    # Backup data from current run for testing
    saveRDS(dataSPIN.ls, file = paste0(export.data, "TMP.RDS"))
  }
}
