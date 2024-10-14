
# Backgrounds
# Activation fraction
{
  # Considering ice as this cutoff size (um) or larger
  cutoff.bins.nm <- bins.nm[which(as.numeric(bins.nm) >= THRESHOLD.Dp.nm[i])]

  # Calculate total particle density larger than cutoff
  tmp.df <- dataALL.df %>%
    mutate(`Total Particles` = rowSums(select(., .dots = all_of(cutoff.bins.nm)))) %>%
    mutate(`Inlet Filter ON` = dataALL.df$`Inlet Filter ON`) %>%
    select(`Inlet Filter ON`, `Total Particles`)

  # Set total backgrounds to NA for when the inlet filter is on
  # Using 0 will sway statistics significantly
  tmp.df <- tmp.df %>%
    mutate(`Total Background` = if_else(`Inlet Filter ON` == 0, `Total Particles`, NA), .after = `Total Particles`)

  # Use a rolling average
  tmp.df <- tmp.df %>%
    mutate(`Total Background` = zoo::rollapplyr(`Total Background`, 1*60, mean, fill = 0))

  # Using a linear interpolation to fill in the time series than averaging gives a higher background ice concentration
  # This gives a more conservative approach in declaring ice counts
  tmp.df <- tmp.df %>%
    mutate(`Total Background` = zoo::na.approx(`Total Background`, na.rm = F))

  # Subtract number density of backgrounds
  tmp.df <- tmp.df %>%
    mutate(`Total Ice` = `Total Particles` - `Total Background`) %>%
    mutate(`Total Ice` = if_else(`Total Ice` > 0, `Total Ice`, 0))

  # Set values when the filter is on to NA
  # Using NA here as it will artifically lower the actual number due to weighting
  tmp.df <- tmp.df %>%
    mutate(`Total Ice` = if_else(`Inlet Filter ON` == 0, 0, `Total Ice`))

  # Apply correction factors
  {
    tmp.df <- tmp.df %>%
      select(!`Inlet Filter ON`)

    tmp.df <- tmp.df*THRESHOLD.CF.nm
    }

  # Calculate activation fraction
  dataAF.df <- cbind(dataALL.df, tmp.df) %>%
    relocate(`Total Particles`, `Total Background`, `Total Ice`, .before = "0.55") %>%
    mutate(`Activation Fraction (%)` = (`Total Ice`/`Total CN`)*100, .after = `Total Ice`) %>%
    mutate(`Lamina Breaks` = as.factor(`Lamina Breaks`))

  ggplot(dataAF.df, aes(x = `Lamina S Ice`, y = `Activation Fraction (%)`, col = `Lamina Breaks`)) +
    geom_point() +
    facet_wrap(~`Class`)

  plot(dataAF.df$`Lamina S Ice`, dataAF.df$`Activation Fraction (%)`, col = as.factor(dataAF.df$`Lamina Breaks`))





}


# ------------------------------------------------------------------------ #
##### SUBSECTION: Activation Analysis #####
#

{
  print(paste0("Plotting Activation Analysis: ", date.c))

  # THRESHOLDS USED
  {
    # WATER ACTIVITY
    {
      THRESHOLD.S.LIQ.nm <- 1.1
    }

    # HOMOGENEOUS FREEZING THRESHOLD
    {
      THRESHOLD.HMG.FRZ.nm <- c(freezing.Dp)
    }

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

    # SIZE LIMIT
    #' These depend on the aerosol testing
    #' Organics will be smaller than mineral dust and should be considered
    {
      THRESHOLD.ALL.nm <- 0.55
      THRESHOLD.ORG.nm <- 2.5
      THRESHOLD.INORG.nm <- 5

      # Plot for each threshold for clarity in choice
      THRESHOLD.Dp.nm <- c(THRESHOLD.ALL.nm, THRESHOLD.ORG.nm, THRESHOLD.INORG.nm)
    }

    # OPC Scattering Limits
    {
      THRESHOLD.SIZE.ALL.nm <- 0
      THRESHOLD.SIZE.ORG.nm <- 3.25
      THRESHOLD.SIZE.INORG.nm <- 3.25

      THRESHOLD.DEPO.ALL.nm <- 0
      THRESHOLD.DEPO.ORG.nm <- 0.3
      THRESHOLD.DEPO.INORG.nm <- 0.4

      THRESHOLD.SIZE.nm <- c(THRESHOLD.SIZE.ALL.nm, THRESHOLD.SIZE.ORG.nm, THRESHOLD.SIZE.INORG.nm)
      THRESHOLD.DEPO.nm <- c(THRESHOLD.DEPO.ALL.nm, THRESHOLD.DEPO.ORG.nm, THRESHOLD.DEPO.INORG.nm)
    }

    # ACTIVATION THRESHOLD
    {
      THRESHOLD.AF.ALL.nm <- 0
      THRESHOLD.AF.ORG.nm <- 0.5
      THRESHOLD.AF.INORG.nm <- 1

      THRESHOLD.AF.nm <- c(THRESHOLD.AF.ALL.nm, THRESHOLD.AF.ORG.nm, THRESHOLD.AF.INORG.nm)
    }
  }

  # Classification
  # This is where the model will acquire the most bias as it is user defined
  {
    # AEROSOL
    {
      # This looks at periods where the chamber walls are set for the same temperature and the inlet filter is open
      class.df <- class.df %>%
        mutate(`Class` = if_else(condition = (class.df$`Warm Wall SP` - class.df$`Cold Wall SP` == 0) & (`Inlet Filter ON` == 1),
                                 true = "Aerosol",
                                 false = `Class`))
    }

    # ICE
    {
      # This looks at periods where the inlet filter is on and 10 um particles are being detected
      # These particles have to be ice from the walls
      class.df <- class.df %>%
        mutate(`Class` = if_else(condition = (`Lamina S Ice` >= `Freezing Threshold`) & (`Large Ice` > 0) & (`Lamina S Liquid` < 1),
                                 true = "Ice",
                                 false = `Class`))
    }

    # DROPLET
    {

      threshold.liq <- quantile(class.df$`Lamina S Liquid`, probs = c(0.99))

      class.df <- class.df %>%
        mutate(`Class` = if_else(condition = (`Lamina S Liquid` >= threshold.liq) & (`Large Ice` > 0),
                                 true = "Droplet",
                                 false = `Class`))
    }
  }

  ONSETS.ls <- list()
  for (i in 1:length(THRESHOLD.Dp.nm)){

    # Backgrounds
    # Activation fraction
    {
      # Considering ice as this cutoff size (um) or larger
      cutoff.bins.nm <- bins.nm[which(as.numeric(bins.nm) >= THRESHOLD.Dp.nm[i])]

      # Calculate total particle density larger than cutoff
      tmp.df <- dataALL.df %>%
        mutate(`Total Particles` = rowSums(select(., .dots = all_of(cutoff.bins.nm)))) %>%
        mutate(`Inlet Filter ON` = dataALL.df$`Inlet Filter ON`) %>%
        select(`Inlet Filter ON`, `Total Particles`)

      # Set total backgrounds to NA for when the inlet filter is on
      # Using 0 will sway statistics significantly
      tmp.df <- tmp.df %>%
        mutate(`Total Background` = if_else(`Inlet Filter ON` == 0, `Total Particles`, NA), .after = `Total Particles`)

      # Using a linear interpolation to fill in the time series than averaging gives a higher background ice concentration
      # This gives a more conservative approach in declaring ice counts
      tmp.df <- tmp.df %>%
        mutate(`Total Background` = zoo::na.approx(`Total Background`, na.rm = F))

      # Use a rolling average
      tmp.df <- tmp.df %>%
        mutate(`Total Background` = zoo::rollapplyr(`Total Background`, 5*60, mean, fill = 0))

      # Subtract number density of backgrounds
      tmp.df <- tmp.df %>%
        mutate(`Total Ice` = `Total Particles` - `Total Background`) %>%
        mutate(`Total Ice` = if_else(`Total Ice` > 0, `Total Ice`, 0))

      # Set values when the filter is on to NA
      # Using NA here as it will artifically lower the actual number due to weighting
      tmp.df <- tmp.df %>%
        mutate(`Total Ice` = if_else(`Inlet Filter ON` == 0, 0, `Total Ice`))

      # Apply correction factors
      {
        tmp.df <- tmp.df %>%
          select(!`Inlet Filter ON`)

        tmp.df <- tmp.df*THRESHOLD.CF.nm
        }

      # Calculate activation fraction
      dataAF.df <- cbind(dataALL.df, tmp.df) %>%
        relocate(`Total Particles`, `Total Background`, `Total Ice`, .before = "0.55") %>%
        mutate(`Activation Fraction (%)` = (`Total Ice`/`Total N`)*100, .after = `Total Ice`)
    }

    # Applying thresholds for filtering
    FINAL.DATA.df <- dataAF.df %>%
      filter(`Lamina S Liquid` < THRESHOLD.S.LIQ.nm) %>%
      filter(`Log Size` > THRESHOLD.SIZE.nm[i]) %>%
      filter(`Depolarization` > THRESHOLD.DEPO.nm[i])



    ONSETS.df <- FINAL.DATA.df %>%
      filter(`Activation Fraction (%)` >= THRESHOLD.AF.nm[i])

    if (nrow(ONSETS.df) != 0) {
      ONSETS.df <- ONSETS.df %>%
        mutate(`Outlier Threshold` = max(min(`Lamina S Ice`, na.rm = T), as.numeric(quantile(`Lamina S Ice`, 0.25)) - (IQR(`Lamina S Ice`)*1.5)), na.rm = T) %>%
        filter(`Lamina S Ice` > `Outlier Threshold`)

      if (nrow(ONSETS.df) != 0) {
        ONSETS.df <- ONSETS.df %>%
          group_by(`Lamina Breaks`, .drop = FALSE) %>%
          summarize(`Counts` = n(), `Onset` = min(`Lamina S Ice`)) %>%
          mutate(`Date` = dates.c[n], .before = everything()) %>%
          mutate(`Threshold Type` = if_else(i == 1, "No Thresholds", NA), .after = `Date`) %>%
          mutate(`Threshold Type` = if_else(i == 2, "Inorganic Aerosol Thresholds", NA), .after = `Date`) %>%
          mutate(`Threshold Type` = if_else(i == 3, "No Thresholds", NA), .after = `Date`)

        ONSETS.ls[[i]] <- ONSETS.df
      }
    }

    if (i == 1){

      WATER.UPTAKE.df <- dataAF.df %>%
        filter(`Lamina S Liquid` < THRESHOLD.S.LIQ.nm) %>%
        filter(`Log Size` > THRESHOLD.SIZE.nm[i]) %>%
        filter(`Activation Fraction (%)` >= 0.5)

      if (nrow(WATER.UPTAKE.df) != 0){

        WATER.UPTAKE.df <- WATER.UPTAKE.df %>%
          mutate(`Outlier Threshold` = max(min(`Lamina S Ice`, na.rm = T), as.numeric(quantile(`Lamina S Ice`, 0.25)) - (IQR(`Lamina S Ice`)*1.5)), na.rm = T) %>%
          filter(`Lamina S Ice` > `Outlier Threshold`)

        if (nrow(WATER.UPTAKE.df) != 0){

          WATER.UPTAKE.df <- WATER.UPTAKE.df %>%
            group_by(`Lamina Breaks`, .drop = FALSE) %>%
            summarize(`Counts` = n(), `Onset` = min(`Lamina S Ice`)) %>%
            mutate(`Date` = dates.c[n], .before = everything()) %>%
            mutate(`Threshold Type` = "Water Uptake", .after = `Date`)
        }
      }
    }

    # Plotting lines for homogeneous freezing at the set diameter
    homogeneous.freezing.thresholds <- FINAL.DATA.df %>%
      filter(`Lamina Breaks` <= -38) %>%
      group_by(`Lamina Breaks`) %>%
      summarise(xleft = min(`Freezing Error Lower`), xright = max(`Freezing Error Upper`))

    caption.df <- as.data.frame(t(data.frame("Minimum Activated Particle Size" = THRESHOLD.Dp.nm[i],
                                             "Activated Fraction Criteria (%)" = THRESHOLD.AF.nm[i],
                                             "Water Saturation Upper Limit (S)" = THRESHOLD.S.LIQ.nm,
                                             "OPC Size Intensity Lower Limit" = THRESHOLD.SIZE.nm[i],
                                             "OPC Depolarization Ratio Lower Limit" = THRESHOLD.DEPO.nm[i],
                                             "Correction Factor" = THRESHOLD.CF.nm,
                                             "Homogeneous Freezing Diameter (um)" = freezing.Dp,
                                             check.names = F)))

    caption.df = caption.df %>%
      mutate(`Label` = paste0(rownames(caption.df)), .before = everything())

    # Plotting
    {
      plot.filename = paste0(export.plot, date.c, "/", date.c, " Activation Analysis", ".pdf")

      pdf(
        file = plot.filename,
        width = 6,
        height = 8,
        bg = "white"
      )

      layout(matrix(c(1, 1, 2, 3), 4, 1, byrow = T))

      # Plot limits
      {
        x.lims = c(1, 1.6)
        x.breaks = seq(1, 1.6, 0.1)

        y.lims = ceiling(range(FINAL.DATA.df$`Activation Fraction (%)`)*10)/10
        y.breaks = pretty(y.lims, 10)

        if (i == 1){
          type = "No Thresholds"

          if (!max(y.lims) >= 0.5){
            y.lims[2] <- 0.5
          }
        }

        if (i == 2){
          type = "Organic Aerosol Thresholds"

          if (!max(y.lims) >= 0.5){
            y.lims[2] <- 0.5
          }
        }

        if (i == 3){
          type = "Inorganic Aerosol Thresholds"

          if (!max(y.lims) >= 1){
            y.lims[2] <- 1
          }
        }
      }

      # Color palette
      {
        group.c <- factor(FINAL.DATA.df$`Lamina Breaks`)

        colors.c <- scales::hue_pal()(nlevels(group.c))

        colors.c <- rev(adjustcolor(colors.c, alpha.f = 0.3))

        if (length(colors.c) > 1){
          palette(colors.c)
        }
      }

      {
        # PLOTTING: Activation Ratio
        {
          par(mar = c(0.5, 5, 4, 5))

          # Generate plot
          plot(FINAL.DATA.df$`Lamina S Ice`,
               FINAL.DATA.df$`Activation Fraction (%)`,
               col = group.c,
               xaxt = "n",
               yaxt = "n",
               xlab = "",
               ylab = "",
               xaxs = 'i',
               cex = 0.8,
               yaxs = 'i',
               xlim = x.lims,
               ylim = y.lims)

          grid(nx = NULL, ny = NULL,
               lty = 1,      # Grid line type
               col = "gray80", # Grid line color
               lwd = 0.5)      # Grid line width

          axis(side = 1, at = x.breaks, labels = x.breaks, tick = T)
          axis(side = 2, at = y.breaks, labels = y.breaks, tick = T)

          Hmisc::minor.tick(nx = 4, ny = 0, tick.ratio = 0.5)

          abline(h = THRESHOLD.AF.nm[i], col = "firebrick4", lwd = 1.5, lty = 2)
          abline(v = ONSETS.df$Onset, col = "black", lwd = 1.5, lty = 1)

          mtext(side = 3, cex = 1, text = paste0(date.c, ": Activation Analysis"), line = 2, outer = F)
          mtext(side = 3, cex = 0.9, text = paste0(type), line = 0.5, outer = F)

          # Legend
          legend("right", inset = c(-0.1, 0.5), xpd = T,
                 legend = levels(group.c),
                 pch = 19, col = colors.c, bty = "n")

          for (i in 1:nrow(homogeneous.freezing.thresholds)){

            rect(xleft = as.numeric(homogeneous.freezing.thresholds[i, "xleft"]), xright = as.numeric(homogeneous.freezing.thresholds[i, "xright"]),
                 ybottom = 0, ytop = 100,
                 col = colors.c[i],
                 border = NA
            )
          }

          mtext(latex2exp::TeX(r'(Activation Fraction (%))'), side = 2, line = 2.5, cex.lab = 0.8, family = "Helvetica")
        }

        # PLOTTING: Depolarization Ratio
        {
          par(mar = c(3, 5, 3, 5))

          y.lims = c(0, 1)
          y.breaks = pretty(y.lims, 4)

          # Generate plot
          plot(dataAF.df$`Lamina S Ice`,
               dataAF.df$`Depolarization`,
               col = group.c,
               xaxt = "n",
               yaxt = "n",
               xlab = "",
               ylab = "",
               xaxs = 'i',
               cex = 0.1,
               yaxs = 'i',
               xlim = x.lims,
               ylim = y.lims)

          grid(nx = NULL, ny = NULL,
               lty = 1,      # Grid line type
               col = "gray80", # Grid line color
               lwd = 0.5)      # Grid line width

          axis(side = 1, at = x.breaks, labels = x.breaks, tick = T, tcl = 0.3)
          axis(side = 2, at = y.breaks, labels = y.breaks, tick = T)

          abline(h = THRESHOLD.DEPO.nm[i], col = "firebrick4", lwd = 1.5, lty = 2)

          Hmisc::minor.tick(nx = 4, ny = 2, tick.ratio = 0.5)

          for (i in 1:nrow(homogeneous.freezing.thresholds)){

            rect(xleft = as.numeric(homogeneous.freezing.thresholds[i, "xleft"]), xright = as.numeric(homogeneous.freezing.thresholds[i, "xright"]),
                 ybottom = -4, ytop = 2,
                 col = colors.c[i],
                 border = NA
            )
          }

          mtext(latex2exp::TeX(r'($Depolarization$)'), side = 2, line = 2.5, cex = 0.8, family = "Helvetica")
          mtext(latex2exp::TeX(r'($S_{ice}$)'), side = 1, line = 2.75, cex = 0.8, family = "Helvetica")
        }

        # TABLE: Parameters
        {
          par(mar = c(5, 4, 5, 2))
          plot.new()

          bg_col <- matrix("white", nrow(caption.df), ncol(caption.df))
          bg_col[c(1, 3, 5, 7), ] <- "grey80"

          plotrix::addtable2plot(0, 0, caption.df,
                                 bty = "o",
                                 bg = bg_col,
                                 xjust = -0.4,
                                 yjust = 0.8,
                                 ypad = 1.05,
                                 display.rownames = F,
                                 display.colnames = F,
                                 hlines = TRUE,
                                 vlines = TRUE)
        }
      }
      }
  }

  dev.off()
}
}
}

{
  # Merge back for exporting
  export.df <- rbindlist(ONSETS.ls) %>%
    filter(`Threshold Type` != "No Thresholds")

  if (nrow(WATER.UPTAKE.df) != 0){
    export.df <- rbind(export.df, WATER.UPTAKE.df)
  }

  # Create a path to export plots with
  export.data.path = paste0(export.data, str_remove_all(dates.c[n], '-'))

  # Check if export path exists
  # If it does not, create it
  if (!dir.exists(export.data.path)) {
    # Create a dated directory to send plots to
    dir.create(export.data.path, mode = "777")
  }

  print(paste0("Exporting Level 1 Data: ", export.data.path))

  # Save data using data.table::fwrite
  data.table::fwrite(export.df, file = paste0(export.data.path, "/", spin.dirs[n], " Summary", ".csv"), showProgress = T)
}
}
}
}



# ------------------------------------------------------------------------ #
##### SUBSECTION: Uncertainty Analysis #####
#

{
  r_um = 0.225

  # Calculate lamina statistics
  tmp.ls <- NULL
  tmp.ls <- lamina.conditions(dataALL.df, dataALL.df$`Lamina Centerline Ratio`, r_um = r_um, t = 60)

  # Determine which indices of the bin vector are above the cutsize
  cutoff.bins.nm <- bins.nm[which(as.numeric(bins.nm) >= 5)]

  tmpx <- dataALL.df %>%
    select(all_of(cutoff.bins.nm)) %>%
    rowSums()

  # Remove averaged background counts
  tmpx.corrected <- tmpx*1.4 - dataALL.df$`Background Ice`*1.4

  # For any rows in which the counts become negative default to 0
  tmpx.corrected[tmpx.corrected < 0] <- 0

  afdfagdgdg <- apply(tmp.ls[[1]][, 2:17], 1, max)

  #
  dataAF.df <- dataALL.df %>%
    mutate(`Total Ice` = tmpx.corrected, .before = "0.55") %>%
    mutate(`Total Ice` = ifelse(`Log S1/P1` >= -0.25, `Total Ice`, 0))

  lamina.df <- cbind("Lamina S Ice" = dataAF.df$`Lamina S Ice`,
                     "Maximum S Ice" = apply(tmp.ls[[1]][, 2:17], 1, max),
                     "Total Ice" = dataAF.df$`Total Ice`,
                     tmp.ls[[1]][, 2:17])

  plot(lamina.df$`Maximum S Ice`, lamina.df$`Total Ice`, pch = 1, cex = 0.1,
       xlab = "Lamina S Ice", ylab = "Total Ice Counts")
  title(main = paste0("Citric Acid Freezing ", date.c))
  points(lamina.df$`Lamina S Ice`, lamina.df$`Total Ice`, pch = 2, col = "red", cex = 0.1)
  abline(v = dataAF.df$`Freezing Threshold`, lwd = 2, col = "blue")
  legend("topleft", legend = c("Maximum S Ice", "Lamina S Ice"), pch = 1:2, col = c("black", "red"), ncol = 1)
}

{
  par(mfrow = c(1,1))

  plot(dataALL.df$`Local Time`, dataALL.df$`Lamina Temp (C)`)

  df <- dataALL.df

  cold.positions = which(str_detect(colnames(df), "C\\d+"))
  warm.positions = which(str_detect(colnames(df), "W\\d+"))

  result1 <- matrix(nrow = nrow(df), ncol = 16)
  result2 <- matrix(nrow = nrow(df), ncol = 16)

  for (r in 1:nrow(df)){

    for (i in 1:16){

      cx = cold.positions[i]
      wx = warm.positions[i]

      cold.x = df[[r, cx]]

      # Thermocouple measurements have the following accuracy
      # ±1C or ±0.75% whichever is greater
      cold.x1 <- 0.75/100
      cold.x2 <- 1 - (cold.x+1)/cold.x

      error = max(c(cold.x1, cold.x2))

      result1[r, i] <- error

      warm.x = df[[r, wx]]

      # Thermocouple measurements have the following accuracy
      # ±1C or ±0.75% whichever is greater
      warm.x1 <- 0.75/100
      warm.x2 <- 1 - (warm.x+1)/warm.x

      error = max(c(warm.x1, warm.x2))

      result2[r, i] <- error
    }
  }

  maximum.cold.error = max(result1, na.rm = T)
  maximum.warm.error = max(result2, na.rm = T)

  average.cold.wall = dataALL.df %>%
    select(all_of(c("C0", "C3", "C6", "C9", "C12"))) %>%
    rowMeans(., na.rm = T)

  average.warm.wall = dataALL.df %>%
    select(all_of(c("W0", "W3", "W6", "W9", "W12"))) %>%
    rowMeans(., na.rm = T)

  dataALL.df$`Cold Wall Mean` = average.cold.wall

  dataALL.df$`Warm Wall Mean` = average.warm.wall



  # This groups by ramp and calculates various statistics on their behavior
  # Using a chi-squared test would be the most accurate
  # Maximum error is the most conservative
  temperature.error <- dataALL.df %>%
    mutate(`Cold Wall Mean` = mean(`C0`:`C12`))
  group_by(`Lamina Breaks`) %>%
    reframe(`Cold Wall Uncertainty 1` = mean(`Cold`)*(0.75/100),
            `Cold Wall Uncertainty 2` = `Lamina Temp (C)` + 1,
            `SD` = sd(`Lamina Temp (C)`),
            `Mean` = mean(`Lamina Temp (C)`),
            `Uncertainty 1` = `Lamina Temp (C)`*(0.75/100),
            `Uncertainty 2` = `Lamina Temp (C)` + 1,
            `SP` = mean(`Lamina Breaks`),
            `Absolute Error` = abs(max(`Lamina Temp (C)` - `Lamina Breaks`)),
            `Relative Error` = (`Absolute Error`/`SP`)*100)

  temperature.error.df <- dataALL.df %>%
    select(`Lamina Breaks`) %>%
    mutate(`Absolute Error` = abs(`Lamina Breaks` - dataALL.df$`Lamina Temp (C)`)) %>%
    mutate(`Lamina Breaks` = as.factor(`Lamina Breaks`))

  # Thermocouple measurements have the following accuracy: ±1C or ±0.75% whichever is greater

  # Convert to long format for easier plotting
  data.long <- temperature.error.df %>%
    pivot_longer(
      cols = !all_of("Lamina Breaks"),
      names_to = "Variable",
      values_to = "Value"
    )

  ggplot(temperature.error.df, aes(x = `Lamina Breaks`, y = `Absolute Error`, group = `Lamina Breaks`)) +
    stat_boxplot(geom = "errorbar", na.rm = T) +
    geom_boxplot()

  ggplot(data = temperature.error.df, aes(x = `Lamina Breaks`, y = `Absolute Error`, group = "Lamina Breaks")) +
    stat_boxplot(geom = "errorbar", na.rm = T)



  cold.positions = which(str_detect(colnames(dataALL.df), "C\\d+"))
  warm.positions = which(str_detect(colnames(dataALL.df), "W\\d+"))

  r_um = 0.300

  tmp.ls <- NULL
  tmp.ls <- lamina.conditions(tmp.df, tmp.df$`Lamina Centerline Ratio`, r_um = r_um, t = 60)
}

















# Loop through list of files and create a dataframe of each
data.ls <- lapply(files.spin, function(x) {
  # Use data.table's fread to read in data
  # Fastest method in reading CSV's and least memory consuming
  # fill MUST equal FALSE with SPIN files
  tmp.df <- data.table::fread(
    paste0(x),
    na.strings = c("", "NA", "NaN"),
    fill = FALSE,
    strip.white = TRUE,
    stringsAsFactors = FALSE
  )

  # Extract date and time from the filename using regular expression
  # (?<=...) is a positive-look behind assertion
  # SPIN files always follow this pattern
  tmp.date = str_extract(x, "(?<=SPIN\\d{3}\\S?)\\d{8}")
  tmp.time = str_extract(x, "(?<=SPIN\\d{3}\\S?\\d{8})\\d{6}")

  # Convert date string in filename to POSIXct POSIXt
  tmp.datetime = as.POSIXct(paste0(tmp.date, " ", tmp.time),
                            format = "%Y%m%d %H%M%S",
                            tz = "America/New_York")

  # Obtain a date for the file at midnight
  tmp.date = as.POSIXct(paste0(tmp.date, " ", "000000"),
                        format = "%Y%m%d %H%M%S",
                        tz = "America/New_York")

  # Create a POSIXct POSIXt object in the dataframe
  tmp.df <- tmp.df %>%
    mutate(`Start Time` = tmp.datetime, .before = everything()) %>%
    mutate(`Local Time` = `Time (sec)` + tmp.date,
           .before = everything())

  return(tmp.df)
})

# Combine data.frames into one large data.frame
# Column numbers must be the same length for each dataframe in order to merge using this function
rawSPIN.df <-
  data.table::rbindlist(data.ls, use.names = T, fill = T)[-1]



























































































+
  scale_color_discrete(labels = c(bquote({0.55 *mu*m <D[p]}<.(group1) *mu*m),
                                  bquote({.(group1) *mu*m <D[p]}<.(group2) *mu*m),
                                  bquote({.(group2) *mu*m <D[p]}< 14.5*mu*m))) +



















































































































































































  # ------------------------------------------------------------------------ #
  ##### SECTION: Time Series #####
#'

{
  bins.indices = which(!is.na(as.numeric(colnames(dataSPIN.df))))

  bins.df <- dataSPIN.df %>%
    select(`Experiment Duration (s)`, `Partition`, all_of(c(bins.indices)))

  bins.long <- reshape2::melt(
    bins.df,
    id.vars = c("Partition", "Experiment Duration (s)"),
    value.name = "Particle Number Concentration",
    variable.name = "Midpoint Diameter"
  )

  title.main <-
    paste0("Spectrometer for Ice Nucleation (SPIN)")
  title.sub <- paste0("Particle Time Series: ", date.c)

  cm.palette = c(
    'white',
    'gray87',
    'gray70',
    'gray48',
    rev(brewer.pal(8, 'YlGnBu')[-1]),
    brewer.pal(8, 'YlOrRd')[-1],
    'red4',
    'black'
  )[-1]

  bins.gg <-
    ggplot(
      bins.long,
      aes(x = `Experiment Duration (s)`, y = `Particle Number Concentration`, color = `Midpoint Diameter`)
    ) +
    geom_line() +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x)
        10 ^ x),
      labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                   .x))
    ) +
    annotation_logticks(outside = F) +
    scale_color_manual(values = cm.palette) +
    coord_cartesian(expand = T, clip = "off") +
    labs(title = title.main,
         subtitle = title.sub) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey80", size = 0.25),
      panel.grid.minor = element_line(colour = "grey80", size = 0.25),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.text.x = element_text(vjust = 0.75, angle = 90),
      axis.title.x = element_text(vjust = -1),
      axis.title.y = element_text(vjust = 2),
      legend.key = element_blank(),
      legend.title = element_text("Midpoint Diameter (um)", hjust = 0.5),
      plot.margin = unit(c(1, 0.1, 1, 0.1), "cm"),
      axis.ticks = element_line(size = 0.5)
    ) + facet_wrap(~ `Partition`, scales = 'free', nrow = 1)

  plot.filename = paste0(export.plot, 'overview/', date.c, '_timeseries', '.png')

  ggsave(
    plot.filename,
    bins.gg,
    width = 16,
    height = 8,
    dpi = resolution.dpi,
    units = "in",
    bg = "#ffffff"
  )
}

# ------------------------------------------------------------------------ #
##### SECTION: Size Distribution Plots #####
#' Custom filled.contour plot
#' Overlayed with supersaturation curve with respect to ice
#' Isothermal conditions annotated in color bars on top axis
#'

{
  cm.palette = c(
    'white',
    'gray90',
    'gray70',
    'gray50',
    rev(brewer.pal(8, 'YlGnBu')[-1]),
    brewer.pal(8, 'YlOrRd')[-1],
    'red4',
    'black'
  )

  substance = "MSSS"

  title = paste0(substance)
  subtitle = paste0("Pre-Cooler Temperature: ", 25, " C")

  time.length = 10
  time.format = "%H:%M"

  date.c <- as.Date(date.c)

  spin.levelplot.png(
    dataSPIN.df,
    cm.palette,
    title,
    subtitle,
    time.type = "Duration",
    time.length,
    time.format,
    date.c,
    export.plot
  )
}

# ---------------------------------------------------------------------- #
##### SECTION: Ice Nuclei Totals #####
#'

ice.totals <- dataSPIN.df %>%
  select(!all_of(as.character(bins.cutoff))) %>%
  select(all_of(which(!is.na(
    as.numeric(colnames(.))
  )))) %>%
  rowSums()

dataSPIN.df <- dataSPIN.df %>%
  mutate(`dN Ice` = ice.totals)

nucleation.df <- dataSPIN.df

nucleation.df$`Lamina S Ice`[which(nucleation.df$`Lamina S Liquid` > 1)] <-
  NA
nucleation.df$`dN Ice`[which(nucleation.df$`Lamina S Liquid` > 1)] <-
  NA

dataSEMS.df <- rawSEMS.df %>%
  select(-c(`Time UTC`)) %>%
  filter(`Local Time` >= sequence.ramps.tm[1])

background_counts.df <- nucleation.df %>%
  filter(`Inlet Filter ON` == 0) %>%
  group_by(`Partition`) %>%
  summarise(across(c(`dN Ice`), ~ mean(., na.rm = T))) %>%
  ceiling(.)

ice.df <- nucleation.df %>%
  select(`Local Time`,
         `Partition`,
         `dN Ice`,
         `Lamina S Ice`,
         `Lamina Temperature`)

for (i in unique(background_counts.df$Partition)) {
  x <- which(ice.df$Partition == background_counts.df$Partition[i])

  ice.df$`dN Ice`[x] <-
    ice.df$`dN Ice`[x] - background_counts.df$`dN Ice`[i]
}



ice.df$`dN Ice`[which(ice.df$`dN Ice` < 0)] <- 0

{
  png(
    paste0(
      export.plot,
      'classification/',
      date.c,
      '_Ice Time Series_',
      '.png'
    ),
    width = 12,
    height = 8,
    units = "in",
    bg = "white",
    res = 400
  )

  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(
    ice.df$`Local Time`,
    ice.df$`dN Ice`,
    type = 'l',
    xlab = "Local Time",
    ylab = "Ice Counts"
  )
  par(new = TRUE)                             # Add new plot
  plot(
    ice.df$`Local Time`,
    ice.df$`Lamina S Ice`,
    pch = 17,
    col = 3,
    type = 'l',
    axes = F,
    xlab = "",
    ylab = ""
  )
  axis(side = 4, at = pretty(range(dataSPIN.df$`Lamina S Ice`)))      # Add second axis
  mtext("SS Ice", side = 4, line = 3)
  title(main = "> 5um Particles")
  abline(v = critical.points, col = "red")
  grid()

  dev.off()
}



ice.df$`Lamina Temperature` <-
  round(ice.df$`Lamina Temperature`, 0)

test <- ice.df %>%
  filter(`Lamina Temperature` == -42)

critical.points = NULL



critical.points <- append(critical.points,
                          test$`Local Time`[last(which(test$`Lamina S Ice` >= 1.469588))])

summary(ice.df)

total.df <- dataSEMS.df %>%
  select(`Local Time`, `dN Total`)

tmp.df <-
  dplyr::full_join(ice.df, total.df, by = "Local Time")

x.full = which(complete.cases(tmp.df))

x.diff = diff(x.full)

window.mean <- NULL
for (i in 1:length(x.diff)) {
  window = (x.full[i] - ceiling(x.diff[i] / 2)):(x.full[i] + ceiling(x.diff[i] /
                                                                       2))
  window = window[window > 0]

  window.mean = append(window.mean,  mean(tmp.df$`dN Ice`[c(window)]) /
                         tmp.df$`dN Total`[x.full[i]])

}

x.full = x.full[-length(x.full)]

tmp.df <-
  data.frame(v1 = tmp.df$`Local Time`[x.full], v2 = window.mean)

setnames(tmp.df,
         old = colnames(tmp.df),
         new = c("Local Time", "Ice Fraction"))

data.df <-
  dplyr::full_join(dataSPIN.df, tmp.df, by = "Local Time")

plot(data.df$`Local Time`, data.df$`dN Ice`)

analysis.df <- data.df %>%
  select(
    `Duration(s)`,
    `Lamina S Ice`,
    `Lamina S Liquid`,
    `Inlet Filter ON`,
    `Lamina Temperature`,
    `dN Ice`,
    `Ice Fraction`
  ) %>%
  filter(`Lamina S Liquid` < 1) %>%
  filter(`Inlet Filter ON` == 1)

# Plot parameters
ice.theme = theme(
  plot.title = element_text(hjust = 0.5, size = 16),
  plot.subtitle = element_text(hjust = 0.5, size = 12),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "grey80", size = 0.25),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA),
  axis.text.x = element_text(vjust = 0.5),
  axis.title.x = element_text(vjust = -1),
  axis.title.y = element_text(vjust = 2),
  plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
  legend.key = element_blank(),
  axis.ticks.length = unit(-1, "mm"),
  axis.ticks = element_line(size = 0.5),
  legend.key.height = unit(2, 'cm'),
  legend.key.width = unit(0.5, 'cm'),
  legend.title.align = 0.5
)

ice.totals.gg <-
  ggplot(data = analysis.df,
         aes(x = `Lamina S Ice`, y = `dN Ice`, colour = `Lamina Temperature`)) +
  geom_point(size = 3, alpha = 0.25) +
  scale_color_viridis(option = "D") +
  labs(title = title,
       subtitle = subtitle,
       colour = "Lamina\nTemperature") +
  ylab(expression("Ice Counts")) + ice.theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0.5, 1, 0.1, 1), "cm")
  )

ice.fraction.gg <-
  ggplot(data = analysis.df,
         aes(x = `Lamina S Ice`, y = `Ice Fraction`, colour = `Lamina Temperature`)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis(option = "D") +
  xlab("Lamina Ice Saturation") +
  ylab(expression("Ice Fraction")) + ice.theme +
  theme(plot.margin = unit(c(0, 1, 0.5, 0.5), "cm"))

ice.gg <-
  ggarrange(
    ice.totals.gg,
    ice.fraction.gg,
    nrow = 2,
    ncol = 1,
    align = "v",
    common.legend = T,
    legend = "right"
  )

plot.filename = paste0(export.plot,
                       'classification/',
                       date.c,
                       '_icenucleation',
                       '.png')

ggsave(
  plot.filename,
  ice.gg,
  width = 8,
  height = 6,
  dpi = resolution.dpi,
  units = "in",
  bg = "#ffffff"
)

# ------------------------------------------------------------------------ #
##### SECTION: Aerosol Generation Totals #####
#'

date.c <- as.Date(unique(as.Date(rawSEMS.df$`Local Time`)))

title = paste0(substance)
subtitle = paste0(date.c)

time.length = 20
time.format = "%H:%M"
time.type = "POSIX"
data = rawSEMS.df
units = "dN"

export.generation = paste0(
  '/home/chrisrapp/Purdue/SPIN/export/plots/generation/',
  date.c,
  " Size Distribution Citric Acid.png"
)

aerosol.levelplot.png(
  data,
  title,
  subtitle,
  time.type,
  time.length,
  time.format,
  date.c,
  units,
  export.generation
)
}
} # CLOSE LOOP
}


# Print duration of loop and print date/status
print(paste0(date.c, ' Elapsed Time (s) ', round((proc.time(
) - pct)[3], 2)))



Temp = seq(-70, -10, length.out = 1000)
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

freezing.line = homogeneous.freezing.solver(seq(-70, -38, length.out = 1000), 0.1, 60)

freezing.error.u = freezing.line + (2.5*p_liquid(Temp)/p_ice(Temp) + 1)/100
freezing.error.l = freezing.line - (2.5*p_liquid(Temp)/p_ice(Temp) + 1)/100

freezing.line[which((freezing.line > line100) == T)] <- NA

test <- data.frame(Temp, SS.ice, SS.liquid, freezing.line, freezing.error.u, freezing.error.l)


data.df <- rbindlist(data.ls, fill = T)

water.df <- data.df %>%
  mutate(`Experiment` = as.character(`Date`)) %>%
  filter(`Threshold Type` == "Water Uptake") %>%
  filter(!`Experiment` %in% c("2022-07-19", "2022-07-21", "2023-09-29", "2023-10-18", "2024-01-24", "2024-03-06"))

water.df$Experiment[water.df$Experiment == "2023-11-02"] <- "MSSS"
water.df$Experiment[water.df$Experiment == "2023-11-03"] <- "MSSS, -70 C"
water.df$Experiment[water.df$Experiment == "2023-11-06"] <- "ESSS"
water.df$Experiment[water.df$Experiment == "2023-11-07"] <- "ESSS, -70 C"
water.df$Experiment[water.df$Experiment == "2023-11-08"] <- "DSSS"
water.df$Experiment[water.df$Experiment == "2024-01-08"] <- "Citric Acid, Thermal"
water.df$Experiment[water.df$Experiment == "2024-01-26"] <- "Citric Acid, Thermal, -30 C"
water.df$Experiment[water.df$Experiment == "2024-02-05"] <- "Citric Acid, Thermal, -70 C"

ggplot(data = test, mapping = aes(x = `Temp`, y = `SS.ice`)) +
  geom_line(data = saturation.lines, aes(x = `Temp`, y = `Value`, group = `Variable`), linetype = 3, color = "gray40", linewidth = 0.5) +
  geom_line(aes(y = line100), color = "black", linewidth = 1) +
  geom_ribbon(data = test, aes(y = `SS.ice`, ymin = freezing.error.l, ymax = freezing.error.u), fill = "grey70", alpha = 0.2) +
  geom_line(aes(y = freezing.line), linetype = 2) +
  geom_ribbon(aes(ymin = line100, ymax = max(`SS.ice`)), fill = "white") +
  scale_y_continuous(breaks = seq(0.8, 1.8, 0.1), minor_breaks = seq(0.8, 1.8, 0.05), expand = c(0, 0), limits = c(1, 1.8),
                     sec.axis = dup_axis()) +
  scale_x_continuous(breaks = seq(-70, -10, 10), expand = c(0, 0), limits = c(-70, -10),
                     sec.axis = dup_axis()) +
  ggprism::annotation_ticks(sides = "trbl", type = "minor", outside = T) +
  labs(title = "Thermodynamic Onset Conditions", subtitle = "Water Uptake") +
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
    panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
    panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.ontop = TRUE,
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
  ) +
  geom_point(data = water.df, aes(y = `Onset`, x = `Lamina Breaks`, color = `Experiment`), shape = 2, size = 2) +

inorganic <- data.df %>%
  mutate(`Experiment` = as.character(`Date`)) %>%
  filter(`Threshold Type` == "Inorganic Aerosol Thresholds")

ggplot(data = test, mapping = aes(x = `Temp`, y = `SS.ice`)) +
  geom_line(data = saturation.lines, aes(x = `Temp`, y = `Value`, group = `Variable`), linetype = 3, color = "gray40", linewidth = 0.5) +
  geom_line(aes(y = line100), color = "black", linewidth = 1) +
  geom_ribbon(data = test, aes(y = `SS.ice`, ymin = freezing.error.l, ymax = freezing.error.u), fill = "grey70", alpha = 0.2) +
  geom_line(aes(y = freezing.line), linetype = 2) +
  geom_ribbon(aes(ymin = line100, ymax = max(`SS.ice`)), fill = "white") +
  scale_y_continuous(breaks = seq(0.8, 1.8, 0.1), minor_breaks = seq(0.8, 1.8, 0.05), expand = c(0, 0), limits = c(1, 1.8),
                     sec.axis = dup_axis()) +
  scale_x_continuous(breaks = seq(-70, -10, 10), expand = c(0, 0), limits = c(-70, -10),
                     sec.axis = dup_axis()) +
  ggprism::annotation_ticks(sides = "trbl", type = "minor", outside = T) +
  geom_point(data = inorganic, aes(y = `Onset`, x = `Lamina Breaks`, color = `Experiment`), size = 2) +
  labs(title = "Thermodynamic Onset Conditions", subtitle = "Inorganic Aerosol Ice Nucleation Thresholds") +
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
    panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
    panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.ontop = TRUE,
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

organic <- data.df %>%
  mutate(`Experiment` = as.character(`Date`)) %>%
  filter(`Threshold Type` == "Organic Aerosol Thresholds") %>%
  filter(!`Experiment` %in% c("2022-07-19", "2022-07-21", "2023-09-29", "2023-10-18", "2024-03-06")) %>%
  filter(`Counts` > 50)

organic$Experiment[organic$Experiment == "2024-02-05"] <- "Citric Acid, Thermal, -70 C"
organic$Experiment[organic$Experiment == "2024-02-07"] <- "Citric Acid, Aqueous, -70 C"


ggplot(data = test, mapping = aes(x = `Temp`, y = `SS.ice`)) +
  geom_line(data = saturation.lines, aes(x = `Temp`, y = `Value`, group = `Variable`), linetype = 3, color = "gray40", linewidth = 0.5) +
  geom_line(aes(y = line100), color = "black", linewidth = 1) +
  geom_ribbon(data = test, aes(y = `SS.ice`, ymin = freezing.error.l, ymax = freezing.error.u), fill = "grey70", alpha = 0.2) +
  geom_line(aes(y = freezing.line), linetype = 2) +
  geom_ribbon(aes(ymin = line100, ymax = max(`SS.ice`)), fill = "white") +
  scale_y_continuous(breaks = seq(0.8, 1.8, 0.1), minor_breaks = seq(0.8, 1.8, 0.05), expand = c(0, 0), limits = c(1, 1.8),
                     sec.axis = dup_axis()) +
  scale_x_continuous(breaks = seq(-70, -10, 10), expand = c(0, 0), limits = c(-70, -10),
                     sec.axis = dup_axis()) +
  ggprism::annotation_ticks(sides = "trbl", type = "minor", outside = T) +
  geom_point(data = organic, aes(y = `Onset`, x = `Lamina Breaks`, color = `Experiment`), size = 2) +
  labs(title = "Thermodynamic Onset Conditions", subtitle = "Organic Aerosol Ice Nucleation Thresholds") +
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
    panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.1),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.1),
    panel.grid.minor = element_line(colour = "grey80", linewidth = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.ontop = TRUE,
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







# ------------------------------------------------------------------------ #
##### SECTION: Partitioning #####
#' Break data down into ramp specific time periods
#' Add factor variable for these ramps
#'

test.df <- dataSPIN.df %>%
  mutate(`Set Point Difference` = round(abs(`Cold Wall SP` - `Warm Wall SP`), 1), .before = everything()) %>%
  mutate(`SP Boolean` = case_when(`Set Point Difference` > 0 ~ 1,
                                  `Set Point Difference` == 0 ~ 0),
         .before = everything()) %>%
  mutate(`Filter ID` = rleid(`Inlet Filter ON`), .before = everything()) %>%
  mutate(`SP ID` = rleid(`SP Boolean`), .before = everything()) %>%
  select(!`SP Boolean`)

df <- data.frame("Local Time" = dataSPIN.df$`Local Time`, check.names = F)

sequence.indexer <- function(df, sequence.name, sequence.label){

  # If the number of starts is only 1
  if (length(sequence.name) == 1){

    df <- df %>%
      mutate("Tmp" = case_when(`Local Time` < sequence.name[1] ~ 0, `Local Time` >= sequence.name[1] ~ 1))

    setnames(df, old = "Tmp", new = sequence.label)
  } else {

    # If there are multiple restarts
    for (t in 2:length(sequence.name)){

      start = sequence.name[t-1]
      end = sequence.name[t]

      df[which(df$`Local Time` >= start & df$`Local Time` < end), sequence.label] <- t - 1

      # This is the last start sequence
      if (t == length(sequence.name)){

        df[which(df$`Local Time` >= sequence.name[t]), sequence.label] <- t
      }
    }

  }

  df[df == 0] <- NA

  return(df)
}

test <- sequence.indexer(df, sequence.start.tm, "Sequence Start ID")
test <- sequence.indexer(test, sequence.icing.tm, "Sequence Icing ID")
test <- sequence.indexer(test, sequence.evrmp.tm, "Sequence Evaporation Ramp ID")
test <- sequence.indexer(test, sequence.ramps.tm, "Sequence Ramp ID")

test.df <- merge(test, test.df, by = "Local Time")

# Partition data into subsets for each experimental ramp
{
  # Determine the end time of ramps (constant SP difference)
  end <- test.df %>%
    filter(`Set Point Difference` > 0) %>%
    group_by(`SP ID`) %>%
    filter(row_number() == n()) %>%
    ungroup() %>%
    mutate("Ramp End Time" = `Local Time`) %>%
    select("Ramp End Time")

  # Determine the start time of ramps (constant SP difference)
  start <- test.df %>%
    filter(`Set Point Difference` > 0) %>%
    group_by(`SP ID`) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate("Ramp Start Time" = `Local Time`) %>%
    select("Ramp Start Time")

  # Create vectors
  end = end[["Ramp End Time"]]
  start = start[["Ramp Start Time"]]

  # Loop to generate the specific time stamps
  c.partitioning.times = NULL
  for (i in 1:length(end)) {
    tmp <- test.df %>%
      filter(`Local Time` >= end[i]) %>%
      filter(`Local Time` <= start[i + 1]) %>%
      filter(`Inlet Filter ON` == 0) # These are the pauses between ramps that measure the backgrounds

    c.partitioning.times = append(c.partitioning.times, tmp$`Local Time`[nrow(tmp) /
                                                                           2])
  }

  # Finalize partitioning vector
  c.partitioning.times = append(first(test.df$`Local Time`), c.partitioning.times)
  c.partitioning.times = append(c.partitioning.times, last(test.df$`Local Time`))

  # Create vectors of times that matches the same length as dataframe
  assignment = NULL
  for (i in 1:(length(c.partitioning.times) - 1)) {
    tmp = which((test.df$`Local Time` >= c.partitioning.times[i]) &
                  (test.df$`Local Time` < c.partitioning.times[i + 1])
    )

    if (i == (length(c.partitioning.times) - 1)) {
      # Correct length of vector
      tmp = append(tmp, nrow(test.df))
    }

    # Ramps that are longer than 15 minutes are kept
    if (length(tmp) > 5 * 60) {
      assignment = append(assignment, rep(i, length(tmp)))
    } else {
      assignment = append(assignment, rep(NA, length(tmp)))
    }
  }

  test.df <- test.df %>%
    mutate(`Partition` = assignment, .before = everything()) %>%
    group_by(`Partition`) %>%
    mutate(`Sequence Duration (s)` = row_number(),
           .before = everything()) %>%
    ungroup()
}



if (length(sequence.ramps.tm) != 0) {

  title.main <-
    paste0("Spectrometer for Ice Nucleation (SPIN)")
  title.sub <- paste0("Automated Partitioning: ", date.c)

  partition.gg <-
    ggplot(test.df, aes(x = `Local Time`)) +
    geom_point(
      aes(y = `SP ID`, color = "SP ID"),
      alpha = 0.1,
      shape = 0,
      size = 0.1
    ) +
    geom_point(
      aes(y = `Filter ID`, color = "Filter ID"),
      alpha = 0.1,
      shape = 0,
      size = 0.1
    ) +
    geom_line(aes(y = `Set Point Difference`, color = "Set Point Absolute Difference"),
              size = 1.5) +
    geom_point(
      aes(y = `Partition`, color = "Partition"),
      alpha = 0.1,
      shape = 0,
      size = 0.1
    ) +
    geom_point(
      aes(y = `Sequence Evaporation Ramp ID`, color = "Evap Ramp"),
      alpha = 0.1,
      shape = 0,
      size = 0.1
    ) +
    geom_point(
      aes(y = `Sequence Start ID`, color = "Sequence Start ID"),
      alpha = 0.1,
      shape = 0,
      size = 0.1
    ) +
    geom_vline(aes(xintercept = sequence.ramps.tm), color = "black") +
    xlab("Local Time") +
    scale_x_datetime(date_breaks = "30 min", date_labels = "%H:%M") +
    scale_y_continuous(n.breaks = 10) +
    labs(title = title.main,
         subtitle = title.sub) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey80", size = 0.25),
      panel.grid.minor = element_line(colour = "grey80", size = 0.25),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_text(vjust = -1),
      axis.title.y = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = "top"
    )

  plot.filename = paste0(export.plot, date.c, '/', 'Partitioning', '.png')

  ggsave(
    plot.filename,
    partition.gg,
    width = 11.5,
    height = 8,
    dpi = resolution.dpi,
    units = "in",
    bg = "#ffffff"
  )
}






# ------------------------------------------------------------------------ #
##### SUBSECTION: Homogeneous Freezing #####
#

{
  # Minimum threshold we would expect to see homogeneous onset
  freezing.Dp <- 0.2

  # Calculate relevant variables
  dataALL.df <- dataALL.df %>%
    mutate(`Freezing Threshold` = homogeneous.freezing.solver(`Lamina Breaks`, freezing.Dp, 60)) %>%
    mutate(`Freezing Error Upper` = `Freezing Threshold` + (2.5*p_liquid(`Lamina Breaks`)/p_ice(`Lamina Breaks`) + 1)/100) %>%
    mutate(`Freezing Error Lower` = `Freezing Threshold` - (2.5*p_liquid(`Lamina Breaks`)/p_ice(`Lamina Breaks`) + 1)/100)
}