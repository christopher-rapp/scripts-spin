##### SCRIPT: SPectrometer for Ice Nucleation (SPIN) Level 2 #####
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
  library(patchwork)

  # Set working directory
  setwd("~/Library/CloudStorage/Box-Box/TAMU Chamber Experiments")
  work.dir <- getwd()

  #' @import
  #' Specify import directory
  import.spin = paste0(work.dir, "/SPIN/export/level1/")
  import.sems = paste0(work.dir, "/SEMS/export/level0/")
  import.pcu = paste0(work.dir, "/PCU/export/level1/")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/SPIN/export/level2/")
  export.plot = paste0(work.dir, "/SPIN/export/plots/")

  #' Function scripts
  source("~/Library/CloudStorage/Box-Box/Organosulfate Proxies/scripts/Functions/functions_aerosols.R")
  source("~/Library/CloudStorage/Box-Box/Organosulfate Proxies/scripts/Functions/functions_microphysics.R")
  source("~/Library/CloudStorage/Box-Box/Organosulfate Proxies/scripts/Functions/functions_spin_plotting.R")
  source("~/Library/CloudStorage/Box-Box/Organosulfate Proxies/scripts/Functions/functions_spin_lamina.R")
  source("~/Library/CloudStorage/Box-Box/Organosulfate Proxies/scripts/Functions/functions_spin_classifier.R")

  # THRESHOLDS USED
  {
    # OPC CORRECTION FACTORS
    #' Uncertainty in counting ice nucleating particles with continuous flow diffusion chambers
    #' https://acp.copernicus.org/articles/17/10855/2017/
    #'
    {
      CF.min = 1.4
      CF.max = 9.5
      CF.mean = 4

      # Select which correction factor to use
      THRESHOLD.CF.nm <- CF.min
    }

    # ACTIVATION THRESHOLD
    {
      THRESHOLD.AF.ALL.nm <- 0
      THRESHOLD.AF.ORG.nm <- 0.5
      THRESHOLD.AF.INORG.nm <- 1

      THRESHOLD.AF.nm <- c(THRESHOLD.AF.ALL.nm, THRESHOLD.AF.ORG.nm, THRESHOLD.AF.INORG.nm)
    }
  }
}

# ---------------------------------------------------------------------------- #
##### SECTION: File Selection #####
#'

{
  {
    # Directories
    spin.files <- list.files(path = import.spin, recursive = FALSE, full.names = FALSE)

    # Pull POSIX style dates from directory names
    spin.dates <- as_date(str_extract(spin.files, '(?<=SPIN_)\\d{8}'), format = '%Y%m%d')
  }

  {
    # Directories
    sems.dirs <- list.dirs(path = import.sems, recursive = FALSE, full.names = FALSE)

    # Pull POSIX style dates from directory names
    sems.dates <- as_date(sems.dirs)
  }

  # Find the intersection of SEMS, PCU, and SPIN data
  dates.c <- as_date(Reduce(intersect, list(spin.dates, sems.dates)))

  # OPTIONAL
  {
    # Directories
    pcu.files <- list.files(path = import.pcu, recursive = FALSE, full.names = T)

    # Pull POSIX style dates from directory names
    pcu.dates <- as_date(str_extract(pcu.files, '(?<=V\\d{1}\\_)\\d{8}'), format = '%Y%m%d')
  }

  # Only proceed with data that has corresponding SEMS data
  {
    # Only keep SEMS dates that match a corresponding SPIN file
    spin.files <- spin.files[which(spin.dates %in% dates.c)]
    sems.dirs <- sems.dirs[which(sems.dates %in% dates.c)]
  }
}

# ---------------------------------------------------------------------------- #
##### SECTION: Level 2 #####
#'

{
  data.export.ls <- NULL
  for (n in 1:length(dates.c)){

    # Select date from loop
    date.c <- dates.c[n]
    print(paste0(date.c, ": Data Aggregation"))

    {
      # Export plot path
      export.plot.path = paste0(export.plot, str_remove_all(date.c, "-"), '/')

      # Obtain all the file paths
      files.spin <- paste0(import.spin, spin.files[n])
      files.sems <- list.files(path = paste0(import.sems, sems.dirs[n]), full.names = T)

      # Check if PCU data is available
      if (date.c %in% pcu.dates){
        files.pcu = pcu.files[which(pcu.dates == date.c)]
      } else {
        files.pcu = NULL
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: SPIN Data #####
      #'

      {
        # Loop through SPIN files and read in
        # SPIN is exported using f.write which preserves time zone in the file
        rawSPIN.df <- fread(files.spin)
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: SEMS Data #####
      #' Read in data from SEMS and extract total counts detected
      #' This currently works for polydisperse settings now
      #' ADD MONODISPERSE OPTION

      {
        data.ls <- lapply(files.sems, function(x){

          # Read in data
          tmp.df <- fread(x) %>%
            mutate(`Local Time` = lubridate::force_tz(`Local Time`, tzone = "America/New_York")) %>%
            mutate(`UTC Time` = lubridate::force_tz(`UTC Time`, tzone = "UTC"))

          return(tmp.df)
        })

        # Merge data
        rawSEMS.df <- rbindlist(data.ls, fill = T)

        tmp.nm <- str_which(colnames(rawSEMS.df), "Total")
        tmp.nm <- append(tmp.nm, str_which(colnames(rawSEMS.df), "Dpg"))
        tmp.nm <- append(tmp.nm, str_which(colnames(rawSEMS.df), "GSD"))

        # Colnames of SEMS data to relocate later
        tmp.c <- colnames(rawSEMS.df)[tmp.nm]

        rawSEMS.df <- rawSEMS.df %>%
          select(`UTC Time`, all_of(tmp.nm))

        # Pad data to 1 Hz
        rawSEMS.df <- padr::pad(rawSEMS.df, interval = "sec", by = "UTC Time")

        # Approximate particle concentrations across all data
        rawSEMS.df <- rawSEMS.df %>%
          mutate(across(.cols = c(2:ncol(.)), ~ zoo::na.approx(.x, na.rm = F)))
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: PCU Data #####
      #'

      if (!is.null(files.pcu)){
        #
        data.ls <- lapply(files.pcu, function(x){

          tmp.df <- fread(x) %>%
            mutate(`Local Time` = lubridate::force_tz(`Local Time`, tzone = "America/New_York"))

          return(tmp.df)
        })

        # Merge and remove the local time string
        rawPCU.df <- rbindlist(data.ls, fill = T) %>%
          select(!`Local Time`)
      } else {
        rawPCU.df <- NULL
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: Aggregate Data #####
      #'
      {
        # Perform joins
        # Right join is focused on only allowing observations from the first join the second if there is a match
        # Inner join is most conservative and only allows observations that match in both
        dataALL.df <- left_join(rawSPIN.df, rawSEMS.df, by = "UTC Time")

        if (!is.null(rawPCU.df)){
          dataALL.df <- left_join(dataALL.df, rawPCU.df, by = "UTC Time")
        }

        # Set lamina breaks to a factor
        dataALL.df <- dataALL.df %>%
          mutate(`Lamina Breaks` = as.factor(`Lamina Breaks`)) %>%
          relocate(all_of(tmp.c), .after = `Depolarization`) %>%
          relocate(c(`Lamina Breaks`, `Depolarization`, `Size Threshold`), .after = `Lamina S Liquid`)

        rm(rawPCU.df, rawSEMS.df, rawSPIN.df, data.ls)
      }
    }

    # ------------------------------------------------------------------------ #
    ##### SECTION: Analysis #####
    #'

    {
      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Activation Analysis #####
      #

      print(paste0(date.c, ": Activation Analysis"))

      # Activation Fraction
      {
        # Subtract number density of backgrounds
        dataALL.df <- dataALL.df %>%
          mutate(`Total Particles` = `Total Particles`*THRESHOLD.CF.nm) %>%
          mutate(`Total Background` = `Total Background`*THRESHOLD.CF.nm) %>%
          mutate(`Total OPC` = `Total Particles` - `Total Background`, .after = `Total Background`) %>%
          mutate(`Total OPC` = if_else(`Total OPC` > 0, `Total OPC`, 0))

        # Set values when the filter is on to NA
        # Using NA here as it will artifically lower the actual number due to weighting
        dataALL.df <- dataALL.df %>%
          mutate(`Total OPC` = if_else(`Inlet Filter ON` == 0, NA, `Total OPC`)) %>%
          mutate(`Activation Fraction (%)` = (`Total OPC`/`Total CN`)*100, .after = `Total OPC`)
      }

      # Activation plots
      {
        plot.filename = paste0(export.plot, date.c, " Activation Analysis.png")
        activation.analysis.png(dataALL.df, plot.filename)
      }
    }

    # ---------------------------------------------------------------------- #
    ##### SUBSECTION: Chamber Analysis #####
    #' Thermocouple pair by pair comparison of expected lamina conditions
    #

    {
      print(paste0(date.c, ": Chamber Analysis"))

      # Smallest diameter that could be a heterogeneous ice nuclei from Vali
      r_um = 0.1

      # Run lamina conditions function on the chamber wall thermocouples.
      # Outputs a list of length 3
      # List 1 - Lamina S ice
      # List 2 - Lamina S liquid
      # List 3 - Homogeneous Freezing Thresholds
      tmp.ls <- NULL
      tmp.ls <- lamina.conditions(dataALL.df, dataALL.df$`Lamina Centerline Ratio`, r_um = r_um, t = 60)

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

      chamber.contour.png(data = tmp.ls[[1]],
                          data2 = tmp.ls[[3]],
                          title = paste0("Spectrometer for Ice Nucleation (SPIN)"),
                          subtitle = paste0("Chamber Ice Saturation: ", date.c),
                          title.key = latex2exp::TeX(r'($S_{ice}$)', bold = T),
                          levels.limits = c(1, 2),
                          col.palette = NULL,
                          n.colors = 20,
                          plot.filename = paste0(export.plot.path, date.c, ' Chamber Ice Saturation', '.png'),
                          plot.width = 12,
                          plot.height = 6,
                          plot.resolution.dpi = 600,
                          caption = caption1.c)

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

      chamber.contour.png(data = tmp.ls[[2]],
                          data2 = NULL,
                          title = paste0("Spectrometer for Ice Nucleation (SPIN)"),
                          subtitle = paste0("Chamber Liquid Saturation: ", date.c),
                          title.key = latex2exp::TeX(r'($S_{liq}$)', bold = T),
                          levels.limits = c(0.5, 1.5),
                          n.colors = 20,
                          col.palette = NULL,
                          plot.filename = paste0(export.plot.path, date.c, ' Chamber Liquid Saturation', '.png'),
                          plot.width = 12,
                          plot.height = 6,
                          plot.resolution.dpi = 600,
                          caption = caption2.c)

      tmp.df <- cbind("UTC Time" = tmp.ls[[1]][, 1], tmp.ls[[1]][, 2:17] - dataALL.df$`Lamina S Ice`)

      caption3.c = paste0("Difference between ice supersaturations derived from
                           individually paired cold-warm wall thermocouples
                           and reported SPIN lamina conditions. Lower supersaturation
                           values at the bottom of the chamber e.g. T12 - T15 are due
                           to the droplet evaporation section.")

      chamber.contour.png(data = tmp.df,
                          data2 = NULL,
                          title = paste0("Spectrometer for Ice Nucleation (SPIN)"),
                          subtitle = paste0("Chamber Observed Ice Saturation Difference: ", date.c),
                          title.key = latex2exp::TeX(r'($\Delta S_{ice}$)', bold = T),
                          levels.limits = c(-0.4, 0.4),
                          col.palette = c(rev(RColorBrewer::brewer.pal(7, 'Reds')), "white", "white", RColorBrewer::brewer.pal(7, 'Blues')),
                          n.colors = NULL,
                          plot.filename = paste0(export.plot.path, date.c, ' Chamber Ice Difference', '.png'),
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
        dataALL.df <- cbind(dataALL.df, tmp.df)
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
        dataALL.df <- cbind(dataALL.df, tmp.df)
      }
    }

    # ------------------------------------------------------------------------ #
    ##### SUBSECTION: Lamina Uncertainty #####
    #

    {
      print(paste0(date.c, ": Lamina Error Analysis"))

      # Calculate error propagation
      error.df <- lamina.error(dataALL.df, position.error = T)

      # Add calculated error propagation values
      dataALL.df <- cbind(dataALL.df, error.df) %>%
        relocate(all_of(colnames(error.df)), .after = `Inlet Filter ON`) %>%
        relocate(c(`Maximum S Ice`, `Maximum S Liq`), .after = `Lamina Temp (C) Error`)

      # Remove error points when there is an issue with the warm wall
      dataALL.df <- dataALL.df %>%
        mutate(`Lamina S Ice Error` = if_else(abs(`Cold Wall SP` - `Warm Wall SP`) == 0 | `Inlet Filter ON` == 0, NA, `Lamina S Ice Error`))

      rm(error.df)
    }

    # ------------------------------------------------------------------------ #
    ##### SECTION: Export Level 2 Data #####

    {
      # Check if export path exists
      # If it does not, create it
      if (!dir.exists(export.data)) {
        # Create a dated directory to send plots to
        dir.create(export.data, mode = "777")
      }

      export.filename = paste0(export.data, "SPIN001_SPIN_", str_remove_all(date.c, '-'), "_level2.csv")

      print(paste0("Exporting Level 2 Data: ", export.filename))

      # Save data using data.table::fwrite
      data.table::fwrite(dataALL.df, file = export.filename, showProgress = T)
    }

    # Export data to a temporary list for optional plotting later
    data.export.ls[[n]] <- dataALL.df
  }

  # Backup data from current run for testing
  saveRDS(data.export.ls, file = paste0(export.data, "TMP.RDS"))
}

  # -------------------------------------------------------------------------- #
  ##### SECTION: Statistics #####
  #'

  {
    # Read in data that has been saved to a list temporarily
    data.export.ls <- readRDS("~/Library/CloudStorage/Box-Box/Organosulfate Proxies/SPIN/export/level2/TMP.RDS")

    summary.ls <- lapply(data.export.ls, function(x){

      print(unique(as.Date(x$`UTC Time`)))

      if (any(str_detect(colnames(x), "PCU"))){

        activated.df <- x %>%
          filter(`Class` == "Ice" | `Class` == "Droplet") %>%
          filter(`Activation Fraction (%)` >= THRESHOLD.AF.nm[2]) %>%
          mutate(`Activation Threshold (%)` = THRESHOLD.AF.nm[2])

        if (nrow(activated.df) != 0){

          # Run statistics of ice activation
          activated.df <- activated.df %>%
            group_by(`Class`, `Lamina Breaks`) %>%
            arrange(., `Lamina S Ice`) %>%
            reframe(`Date` = unique(as.Date(`UTC Time`)),
                    `Class` = unique(`Class`),
                    `Size Threshold` = unique(`Size Threshold`),
                    `Maximum Activation (%)` = max(`Activation Fraction (%)`, na.rm = T),
                    `Onset RH` = first(`Lamina S Ice`),
                    `Onset RH Error` = first(`Lamina S Ice Error`),
                    `Onset Temperature` = first(`Lamina Temp (C)`),
                    `Onset Temperature Error` = first(`Lamina Temp (C) Error`),
                    `Dpg (um)` = mean(`Dpg`, na.rm = T),
                    `Geometric Standard Deviation` = mean(`GSD`, na.rm = T),
                    `CPC` = mean(`Total CN`, na.rm = T),
                    `PCU Lamina` = mean(`TC1`, na.rm = T),
                    `PCU Lamina Error` = max(`PCU Lamina Error`, na.rm = T),
                    `PCU RH` = mean(`RH2`, na.rm = T),
                    `PCU RH Error` = max(`RH2 Error`, na.rm = T),
                    `Inlet RH` = mean(`RH1`, na.rm = T),
                    `Inlet RH Error` = max(`RH1 Error`, na.rm = T))

          unactivated.df <- NULL

        } else {

          activated.df <- NULL

          unactivated.df <- x %>%
            group_by(`Lamina Breaks`) %>%
            reframe(`Date` = unique(as.Date(`UTC Time`)),
                    `Class` = NA,
                    `Size Threshold` = unique(`Size Threshold`),
                    `Maximum Activation (%)` = max(`Activation Fraction (%)`, na.rm = T),
                    `Onset RH` = NA,
                    `Onset RH Error` = NA,
                    `Onset Temperature` = NA,
                    `Onset Temperature Error` = NA,
                    `Dpg (um)` = mean(`Dpg`, na.rm = T),
                    `Geometric Standard Deviation` = mean(`GSD`, na.rm = T),
                    `CPC` = mean(`Total CN`, na.rm = T),
                    `PCU Lamina` = mean(`TC1`, na.rm = T),
                    `PCU Lamina Error` = max(`PCU Lamina Error`, na.rm = T),
                    `PCU RH` = mean(`RH2`, na.rm = T),
                    `PCU RH Error` = max(`RH2 Error`, na.rm = T),
                    `Inlet RH` = mean(`RH1`, na.rm = T),
                    `Inlet RH Error` = max(`RH1 Error`, na.rm = T))
        }

        export.df <- rbind(activated.df, unactivated.df)
      }
    })

    export.df <- rbindlist(summary.ls, use.names = T, fill = T)

    write.csv(export.df, file = paste0(export.data, "summary_level2.csv"))
  }

  # -------------------------------------------------------------------------- #
  ##### SECTION: Merge Literature and Experiment List #####
  #'

  {
    compounds.df <- readxl::read_xlsx('~/Library/CloudStorage/Box-Box/Organosulfate Proxies/README_Compounds.xlsx', na = "NA")
    literature.df <- readxl::read_xlsx('~/Library/CloudStorage/Box-Box/Organosulfate Proxies/README_Literature.xlsx', na = "NA")

    data.df <- right_join(export.df, compounds.df, by = "Date") %>%
      select(`Date`, `Compound`, `Tg,Dry (K)`, `Heating Temperature`, everything()) %>%
      arrange(`Date`) %>%
      filter(!is.na(`Lamina Breaks`))

    data.df[sapply(data.df, is.infinite)] <- NA
    data.df[sapply(data.df, is.nan)] <- NA

    for (i in 1:nrow(data.df)){

      if (!is.na(data.df$`Tg,Dry (K)`[i]) & !is.na(data.df$`PCU RH`[i])){

        tmp <- Tg.calc(data.df$`Dpg (um)`[i]*1000, data.df$`PCU RH`[i], k = data.df$Kappa[i], kGT = 2.5, TgDry = data.df$`Tg,Dry (K)`[i], Error = data.df$Error[i])

        data.df$Tg[i] = tmp[1]
        data.df$TgError[i] = tmp[2]
      } else {
        data.df$Tg[i] = NA
        data.df$TgError[i] = NA
      }
    }

    data.df <- data.df %>%
      select(`Date`, `Compound`, `Tg,Dry (K)`, `Error`, `Tg`, `TgError`, `Heating Temperature`, everything())

    other.df <- data.df %>%
      filter(`Class` == "Droplet" | is.na(`Class`)) %>%
      mutate(`Phase` = "Liquid")

    data.df <- data.df %>%
      filter(`Class` == "Ice") %>%
      mutate(`Phase` = "Glassy")

    data.df <- rbind(data.df, literature.df, use.names = T, fill = T)

    data.df <- data.df %>%
      mutate(`Phase` = if_else(is.na(`Phase`), "Glassy", `Phase`)) %>%
      mutate(`Class` = if_else(is.na(`Class`), "Ice", `Class`))

    data.df <- rbind(data.df, other.df, fill = T) %>%
      filter(!str_detect(`Compound`, "ammonium")) %>%
      mutate(`Class` = factor(`Class`, levels = c("Ice", "Droplet"))) %>%
      mutate(`Phase` = factor(`Phase`, levels = c("Liquid", "Glassy")))

    data.df <- data.df %>%
      mutate(`Label` = if_else(is.na(`Source`), `Compound`, paste0(`Compound`, " (", `Source`, ")")))

    # This removes any of the data points from the table that correspond to literature
    export.df <- data.df %>%
      filter(!is.na(`Date`)) %>%
      arrange(`Class`, `Compound`) %>%
      mutate(`Generation Method` = if_else(!is.na(`Heating Temperature`), paste0(`Heating Temperature`, " (°C), ", Flowrate, " L min-1"), "Atomizer")) %>%
      mutate(`Tg,Dry (°C)` = paste0(round(`Tg,Dry (K)` - 273.15, 0), " ± ", `Error`), ) %>%
      mutate(`Tg (°C)` = paste0(round(`Tg` - 273.15, 0), " ± ", round(`TgError`, 0))) %>%
      mutate(across(`Onset RH`:`Onset Temperature Error`, \(x) round(x, 2))) %>%
      mutate(`Onset Sice` = paste0(`Onset RH`, " ± ", `Onset RH Error`)) %>%
      mutate(`Onset T (°C)` = paste0(`Onset Temperature`, " ± ", `Onset Temperature Error`)) %>%
      mutate(`PCU T (°C)` = paste0(round(`PCU Lamina`, 2), " ± ", round(`PCU Lamina Error`, 2))) %>%
      mutate(`PCU RH (%)` = paste0(round(`PCU RH`, 0), " ± ", round(`PCU RH Error`, 1))) %>%
      mutate(`PCU Inlet RH (%)` = paste0(round(`Inlet RH`, 0), " ± ", round(`Inlet RH Error`, 1))) %>%
      mutate(`Dpg (um)` = round(`Dpg (um)`, 3)) %>%
      mutate(`σg` = round(`Geometric Standard Deviation`, 2)) %>%
      mutate(`CPC (n cm^-3)` = round(CPC, 0)) %>%
      relocate(`Generation Method`, .after = `Compound`) %>%
      select(c(`Date`, Compound, `Generation Method`, `Tg,Dry (°C)`, `Tg (°C)`,
               `PCU T (°C)`, Class, `Onset Sice`, `Onset T (°C)`,
               `PCU RH (%)`, `PCU Inlet RH (%)`, `Dpg (um)`, `σg`, `CPC (n cm^-3)`))

    write.csv(export.df, file = paste0(export.data, "Table.csv"))

  }


  # -------------------------------------------------------------------------- #
  ##### SECTION: Thermodynamic Onset Figure #####
  #'

  # -------------------------------------------------------------------------- #
  ##### SECTION: Light Scattering Figures #####
  #'
  #'







  all.df <- rbindlist(data.export.ls, fill = T)

  {
    All.df <- all.df %>%
      filter(`Inlet Filter ON` == 1) %>%
      mutate(`Date` = as.Date(`UTC Time`))

    All.df <- right_join(All.df, compounds.df, by = "Date")

    All.df <- All.df %>%
      filter(str_detect(`Compound`, "bisulfate")) %>%
      filter(str_detect(`Compound`, "ammonium")) %>%
      mutate(`Saturated` = if_else(`Lamina S Liquid` >= 1, "Saturated", "Subsaturated")) %>%
      mutate(`Lamina Breaks` = factor(`Lamina Breaks`, levels = c(-35, -40, -45)))
  }

  gg1 <- ggplot(All.df, aes(y = `Depolarization`, x = after_stat(scaled), fill = `Lamina Breaks`)) +
    geom_density(alpha = 0.4) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
    ylab(label = latex2exp::TeX(r'($\delta_{SPIN}$)')) +
    xlab(label = "Normalized Density") +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
          panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
          panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
          axis.title.x = element_text(vjust = -1),
          axis.title.y = element_text(vjust = 4),
          axis.ticks.length = unit(2, "mm"),
          axis.minor.ticks.length = unit(1, "mm"),
          legend.title = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  gg2 <- ggplot(All.df, aes(x = `Lamina S Ice`, y = `Depolarization`, color = `Lamina Breaks`)) +
    geom_point(size = 0.01, alpha = 0.1) +
    geom_smooth(alpha = 0.5) +
    scale_x_continuous(breaks = seq(1, 1.6, 0.1), limits = c(1, 1.6), expand = c(0.001, 0.001)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0.001, 0.001)) +
    xlab(label = latex2exp::TeX(r'($S_{ice}$)')) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
          panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
          panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_blank(),
          axis.ticks.length = unit(2, "mm"),
          axis.minor.ticks.length = unit(1, "mm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
    guides(color = FALSE)

  combined <- (gg1 & theme(plot.tag.position  = c(.935, .96))) -
    (gg2 & theme(plot.tag.position  = c(.785, .96))) +
    plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect") & theme(legend.position = "right")

  plot.filename = paste0(export.plot, "Light.png")

  ggsave(plot.filename, combined, width = 10, height = 6)
}




{
  Temp = seq(-75, -10, length.out = 1000)
  SS.ice = seq(0.80, 1.80, length.out = 1000)
  SS.liquid = SS.ice*(p_ice(Temp)/p_liquid(Temp))

  line70 = (70*p_liquid(Temp)/p_ice(Temp))/100
  line80 = (80*p_liquid(Temp)/p_ice(Temp))/100
  line90 = (90*p_liquid(Temp)/p_ice(Temp))/100
  line100 = (100*p_liquid(Temp)/p_ice(Temp))/100

  saturation.lines = data.frame(SS.ice, Temp, line70, line80, line90)

  # Convert to long format for easier plotting
  saturation.lines <- saturation.lines %>%
    pivot_longer(
      cols = !c("SS.ice", "Temp"),
      names_to = "Variable",
      values_to = "Value"
    )

  freezing.line = homogeneous.freezing.solver(seq(-75, -10, length.out = 1000), 0.2, 60)

  freezing.error.u = freezing.line + (2.5*p_liquid(Temp)/p_ice(Temp) + 1)/100
  freezing.error.l = freezing.line - (2.5*p_liquid(Temp)/p_ice(Temp) + 1)/100

  test <- data.frame(Temp, SS.ice, SS.liquid, line100, freezing.line, freezing.error.u, freezing.error.l)

  back.gg <- ggplot(data = test, mapping = aes(x = `Temp`, y = `SS.ice`)) +
    geom_ribbon(aes(y = `SS.ice`, ymin = freezing.error.l, ymax = freezing.error.u), fill = "grey70", alpha = 0.2) +
    geom_line(aes(y = freezing.line), linetype = 2) +
    geom_ribbon(aes(ymin = line100, ymax = max(`SS.ice`)), fill = "white") +
    geom_vline(xintercept = seq(-75, -10, by = 10), col = "gray80", linewidth = 0.1) +
    geom_hline(yintercept = seq(1, 1.8, by = 0.1), col = "gray80", linewidth = 0.1) +
    geom_vline(xintercept = seq(-75, -10, by = 5), col = "gray90", linewidth = 0.1) +
    geom_hline(yintercept = seq(1, 1.8, by = 0.05), col = "gray90", linewidth = 0.1) +
    geom_line(data = saturation.lines, aes(x = `Temp`, y = `Value`, group = `Variable`), linetype = 3, color = "gray40", linewidth = 0.5) +
    geom_line(aes(y = line100), color = "black", linewidth = 1) +
    scale_y_continuous(breaks = seq(0.8, 1.8, 0.1), minor_breaks = seq(0.8, 1.8, 0.05), expand = c(0, 0), limits = c(1, 1.8),
                       sec.axis = dup_axis()) +
    scale_x_continuous(breaks = seq(-75, -10, 10), expand = c(0, 0), limits = c(-75, -10)) +
    ggprism::annotation_ticks(sides = "rbl", type = "minor", outside = T) +
    xlab("Temperature (\u00B0C)") +
    ylab(bquote(S[ice])) +
    coord_cartesian(clip = "off") +
    annotate(geom = "label", x = -42, y = 1.05, label = as.character(expression(S[liq] == 0.7)), parse = T, fill = "white", label.size = NA, size = 10/.pt, color = "gray40") +
    annotate(geom = "label", x = -28, y = 1.05, label = as.character(expression(S[liq] == 0.8)), parse = T, fill = "white", label.size = NA, size = 10/.pt, color = "gray40") +
    annotate(geom = "label", x = -15, y = 1.05, label = as.character(expression(S[liq] == 0.9)), parse = T, fill = "white", label.size = NA, size = 10/.pt, color = "gray40") +
    theme(
      plot.title = element_text(),
      plot.subtitle = element_text(color = "gray25"),
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.title.x.top = element_blank(),
      axis.text.x.top = element_blank(),
      axis.title.x = element_text(),
      axis.text.x = element_text(),
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.ticks.y = element_line(linewidth = 0.5),
      axis.ticks.length = unit(2, "mm"),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
}


hash.width.x = diff(layer_scales(back.gg)$x$range$range) / 40
hash.width.y = diff(layer_scales(back.gg)$y$range$range) / 40

# Colorblind accessible color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

back.gg + geom_errorbar(
  data = data.df,
  aes(
    x = `Onset Temperature`,
    y = `Onset RH`,
    ymin = `Onset RH` - `Onset RH Error`,
    ymax = `Onset RH` + `Onset RH Error`,
    color = `Label`
  ),
  width = hash.width.x,
  linewidth = 0.25
) +
  geom_errorbarh(
    data = data.df,
    aes(
      x = `Onset Temperature`,
      y = `Onset RH`,
      xmin = `Onset Temperature` - `Onset Temperature Error`,
      xmax = `Onset Temperature` + `Onset Temperature Error`,
      color = `Label`
    ),
    height = hash.width.y,
    linewidth = 0.25
  ) +
  geom_point(
    data = data.df,
    aes(
      y = `Onset RH`,
      x = `Onset Temperature`,
      color = `Label`,
      shape = `Phase`
    ),
    alpha = 0.75,
    size = 3
  ) +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio = 1, panel.spacing = unit(10, 'mm')) +
  facet_wrap(~`Class`)

plot.filename = paste0(export.plot, "Summary.png")
ggsave(plot.filename, width = 12, height = 8)



library(insight)

export_table(format_table(export.df))




