UFAx_workflow <- function(spreadsheet) {
  ##
  gc()
  closeAllConnections()
  ##
  exhaustive_chemical_enumeration_annotated_table <- NULL
  ##
  ##############################################################################
  ##
  initiation_time <- Sys.time()
  message("Initiated testing the spreadsheet consistency!")
  ##
  checkpoint_parameter <- FALSE
  #
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_exECS <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      message("The UFAx spreadsheet was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        spreadsheet <- readxl::read_xlsx(spreadsheet, sheet = "exhaustive_chemical_enumeration")
        PARAM_exECS <- cbind(spreadsheet[, 2], spreadsheet[, 4])
        checkpoint_parameter <- TRUE
      } else {
        message("The UFAx spreadsheet not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      message("The UFAx spreadsheet was not produced properly!")
    }
  } else {
    message("The UFAx spreadsheet was not produced properly!")
  }
  ##
  ##############################################################################
  ##
  if (checkpoint_parameter) {
    ##
    x0001 <- which(PARAM_exECS[, 1] == "exECS0001")
    inputHRMSfolderPath <- gsub("\\", "/", PARAM_exECS[x0001, 2], fixed = TRUE)
    PARAM_exECS[x0001, 2] <- inputHRMSfolderPath
    if (dir.exists(inputHRMSfolderPath)) {
      MSfileName <- PARAM_exECS[which(PARAM_exECS[, 1] == "exECS0002"), 2]
      if (!file.exists(paste0(inputHRMSfolderPath, "/", MSfileName))) {
        checkpoint_parameter <- FALSE
        message("ERROR!!! Problem with exECS0002! HRMS is not available!")
      }
    } else {
      checkpoint_parameter <- FALSE
      message("ERROR!!! Problem with exECS0001! Folder of HRMS file is not available!")
    }
    ##
    x0003 <- which(PARAM_exECS[, 1] == "exECS0003")
    addressPeaklist <- gsub("\\", "/", PARAM_exECS[x0003, 2], fixed = TRUE)
    PARAM_exECS[x0003, 2] <- addressPeaklist
    if (dir.exists(addressPeaklist)) {
      peaklistFileName <- PARAM_exECS[which(PARAM_exECS[, 1] == "exECS0004"), 2]
      addressPeaklistFileName <- paste0(addressPeaklist, "/", peaklistFileName)
      if (!file.exists(addressPeaklistFileName)) {
        checkpoint_parameter <- FALSE
        message("ERROR!!! Problem with exECS0004! peaklist is not available!")
      }
    } else {
      checkpoint_parameter <- FALSE
      message("ERROR!!! Problem with exECS0003! Folder of peaklist is not available!")
    }
    peaklist <- IDSL.IPA::loadRdata(addressPeaklistFileName)
    noPeaks <- dim(peaklist)[1]
    ##
    exECS0005 <- PARAM_exECS[which(PARAM_exECS[, 1] == "exECS0005"), 2]
    if (is.na(exECS0005)) {
      checkpoint_parameter <- FALSE
      message("ERROR!!! Problem with exECS0005! This parameter should be `All` or a vector of indices!")
    } else {
      if (gsub(" ", "", tolower(exECS0005)) == "all") {
        selectedIDpeaklist <- 1:noPeaks
        message("The enitre 12C m/z values in the peaklist were placed in the processing row! Annotated molecular formulas for peak IDs are kept in the 'log_exECS_annotation_' folder!")
      } else {
        selectedIDpeaklist <- tryCatch(eval(parse(text = paste0("c(", exECS0005, ")"))), error = function(e) {NULL})
        if (is.null(selectedIDpeaklist) | (max(selectedIDpeaklist) > noPeaks)) {
          checkpoint_parameter <- FALSE
          message("ERROR!!! Problem with exECS0005! The range of indices are out of the peaklist dimension!")
        } else {
          message("The following peak IDs were selected for processing: ")
          for (id in 1:length(selectedIDpeaklist)) {
            message(paste0(selectedIDpeaklist[id], " - ", peaklist[selectedIDpeaklist[id], 3],  " - ", peaklist[selectedIDpeaklist[id], 8]))
          }
        }
      }
    }
    ##
    x0006 <- which(PARAM_exECS[, 1] == "exECS0006")
    output_path <- gsub("\\", "/", PARAM_exECS[x0006, 2], fixed = TRUE)
    PARAM_exECS[x0006, 2] <- output_path
    if (!dir.exists(output_path)) {
      tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w) {warning("Problem with exECS0006! R cannot create the folder!")})
      if (!dir.exists(output_path)) {
        checkpoint_parameter <- FALSE
      }
    }
    ##
    number_processing_threads <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0007'), 2])
    if (length(number_processing_threads) == 0) {
      message("ERROR!!! Problem with exECS0007! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          message("ERROR!!! Problem with exECS0007! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        message("ERROR!!! Problem with exECS0007! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    maxR13C <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0008'), 2])
    if (length(maxR13C) == 0) {
      message("ERROR!!! Problem with exECS0008! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (maxR13C <= 0) {
        message("ERROR!!! Problem with exECS0008! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    ########################### Chemical space #################################
    ############################################################################
    B_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0009'), 2])
    if (length(B_MAX) == 0) {
      message("ERROR!!! Problem with exECS0009! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (B_MAX < 0) {
        message("ERROR!!! Problem with exECS0009! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Br_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0010'), 2])
    if (length(Br_MAX) == 0) {
      message("ERROR!!! Problem with exECS0010! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Br_MAX < 0) {
        message("ERROR!!! Problem with exECS0010! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Cl_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0011'), 2])
    if (length(Cl_MAX) == 0) {
      message("ERROR!!! Problem with exECS0011! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Cl_MAX < 0) {
        message("ERROR!!! Problem with exECS0011! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    S_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0012'), 2])
    if (length(S_MAX) == 0) {
      message("ERROR!!! Problem with exECS0012! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (S_MAX < 0) {
        message("ERROR!!! Problem with exECS0012! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Si_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0013'), 2])
    if (length(Si_MAX) == 0) {
      message("ERROR!!! Problem with exECS0013! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Si_MAX < 0) {
        message("ERROR!!! Problem with exECS0013! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    N_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0014'), 2])
    if (length(N_MAX) == 0) {
      message("ERROR!!! Problem with exECS0014! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (N_MAX < 0) {
        message("ERROR!!! Problem with exECS0014! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    As_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0015'), 2])
    if (length(As_MAX) == 0) {
      message("ERROR!!! Problem with exECS0015! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (As_MAX < 0) {
        message("ERROR!!! Problem with exECS0015! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    I_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0016'), 2])
    if (length(I_MAX) == 0) {
      message("ERROR!!! Problem with exECS0016! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (I_MAX < 0) {
        message("ERROR!!! Problem with exECS0016! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    O_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0017'), 2])
    if (length(O_MAX) == 0) {
      message("ERROR!!! Problem with exECS0017! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (O_MAX < 0) {
        message("ERROR!!! Problem with exECS0017! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    P_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0018'), 2])
    if (length(P_MAX) == 0) {
      message("ERROR!!! Problem with exECS0018! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (P_MAX < 0) {
        message("ERROR!!! Problem with exECS0018! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Na_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0019'), 2])
    if (length(Na_MAX) == 0) {
      message("ERROR!!! Problem with exECS0019! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (Na_MAX < 0) {
        message("ERROR!!! Problem with exECS0019! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    K_MAX <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0020'), 2])
    if (length(K_MAX) == 0) {
      message("ERROR!!! Problem with exECS0020! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (K_MAX < 0) {
        message("ERROR!!! Problem with exECS0020! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    ipw <- PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0021'), 2]
    if (ipw == "[M+H/K/Na]" | ipw == "[M+H]" | ipw == "[M+Na]" | ipw == "[M+K]") {
      ipw_n <- -1
    } else if (ipw == "[M-H]") {
      ipw_n <- +1
    } else if (ipw == "[M]") {
      ipw_n <- 0
    } else {
      message("ERROR!!! Problem with exECS0021! This parameter should be any of '[M+H/K/Na]', '[M-H]', '[M]'!")
      checkpoint_parameter <- FALSE
    }
    ##
    peak_spacing <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0022'), 2])
    if (length(peak_spacing) == 0) {
      message("ERROR!!! Problem with exECS0022! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (peak_spacing < 0) {
        message("ERROR!!! Problem with exECS0022! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    intensity_cutoff_str <- PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0023'), 2]
    if (length(intensity_cutoff_str) == 0) {
      message("ERROR!!! Problem with exECS0023!")
      checkpoint_parameter <- FALSE
      ##
      c <- 5
      b <- 5
      br <- 5
      cl <- 5
      k <- 5
      s <- 5
      se <- 5
      si <- 5
      ##
      checkStrInt <- FALSE
      tryCatch(eval(parse(text = intensity_cutoff_str)), error = function(e) {checkStrInt <- TRUE})
      if (checkStrInt) {
        message("ERROR!!! Problem with exECS0023!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    UFA_IP_memeory_variables <- tryCatch(eval(parse(text = paste0("c(", PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0024'), 2], ")"))), error = function(e) {NULL})
    if (length(UFA_IP_memeory_variables) != 3) {
      message("ERROR!!! Problem with exECS0024! This parameter should be a vector of three positive numbers!")
      checkpoint_parameter <- FALSE
    }
    ##
    maxHalogenCounts <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0025'), 2])
    if (length(maxHalogenCounts) == 0) {
      message("ERROR!!! Problem with exECS0025! This parameter should be a number!")
      checkpoint_parameter <- FALSE
    } else {
      if (maxHalogenCounts < 0) {
        message("ERROR!!! Problem with exECS0025! This parameter should be a number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    Na_K_x <- PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0026'), 2]
    if (Na_K_x == "1" | tolower(Na_K_x) == "t" | tolower(Na_K_x) == "true") {
      NaKruleCheck <- TRUE
    } else if (Na_K_x == "" | Na_K_x == " " | Na_K_x == "0" | tolower(Na_K_x) == "f" | tolower(Na_K_x) == "false") {
      NaKruleCheck <- FALSE
    } else {
      message("ERROR!!! Problem with exECS0026!")
      checkpoint_parameter <- FALSE
    }
    ##
    MaxR13C <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0027'), 2]) + maxR13C
    if (length(MaxR13C) == 0) {
      message("ERROR!!! Problem with exECS0027! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (MaxR13C <= 0) {
        message("ERROR!!! Problem with exECS0027! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    ################# Molecular formula annotation criteria ####################
    ############################################################################
    massAccuracy <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0028'), 2])
    if (length(massAccuracy) == 0) {
      message("ERROR!!! Problem with exECS0028! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (massAccuracy > 0.01) {
        message("ERROR!!! Problem with exECS0028! Mass accuracy must be below `0.01 Da`")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    maxNEME <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0029'), 2])
    if (length(maxNEME) == 0) {
      message("ERROR!!! Problem with exECS0029! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (maxNEME <= 0) {
        message("ERROR!!! Problem with exECS0029! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minPCS <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0030'), 2])
    if (length(minPCS) == 0) {
      message("ERROR!!! Problem with exECS0030! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (minPCS <= 0) {
        message("ERROR!!! Problem with exECS0030! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minNDCS <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0031'), 2])
    if (length(minNDCS) == 0) {
      message("ERROR!!! Problem with exECS0031! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (minNDCS <= 0) {
        message("ERROR!!! Problem with exECS0031! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minRCS <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0032'), 2])
    if (length(minRCS) == 0) {
      message("ERROR!!! Problem with exECS0032! This parameter should be between 0-100!")
      checkpoint_parameter <- FALSE
    } else {
      if ((minRCS < 0) | (minRCS > 100)) {
        message("ERROR!!! Problem with exECS0032! This parameter should be between 0-100!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    scoreCoefficients <- tryCatch(eval(parse(text = PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0033'), 2])), error = function(e) {NULL})
    if (is.null(scoreCoefficients)) {
      message("ERROR!!! Problem with exECS0033! This parameter should be a vector of five positive numbers!")
      checkpoint_parameter <- FALSE
    } else {
      if (length(scoreCoefficients) != 5) {
        message("ERROR!!! Problem with exECS0033! This parameter should be a vector of five positive numbers!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    maxAllowedNumberHits <- as.numeric(PARAM_exECS[which(PARAM_exECS[, 1] == "exECS0034"), 2])
    if (length(maxAllowedNumberHits) == 0) {
      message("ERROR!!! Problem with exECS0034! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (maxAllowedNumberHits > 0) {
        dev.offCheck <- TRUE
        while (dev.offCheck) {
          dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
        }
        ##
        exportSpectraCheck <- TRUE
        outputProfileSpectra <- paste0(output_path, "/UFAx_spectra")
        if (!dir.exists(outputProfileSpectra)) {
          dir.create(outputProfileSpectra, recursive = TRUE)
        }
        message("UFAx_spectra comparison plots with theoretical isotopic profiles are stored in the `UFAx_spectra` folder!")
      } else {
        exportSpectraCheck <- FALSE
        exportSpectraParameters <- NULL
      }
    }
    ##
    exECS0035 <- PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0035'), 2]
    if (length(exECS0035) == 0) {
      message("ERROR!!! Problem with exECS0035!")
      checkpoint_parameter <- FALSE
    } else {
      if (!(tolower(exECS0035) == "yes" | tolower(exECS0035) == "no")) {
        message("ERROR!!! Problem with exECS0035!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (tolower(exECS0035) == "yes") {
      IonPathways <- tryCatch(eval(parse(text = paste0("c(", PARAM_exECS[which(PARAM_exECS[, 1] == 'exECS0036'), 2], ")"))), error = function(e) {NULL})
      if (is.null(IonPathways)) {
        message("ERROR!!! Problem with exECS0036!")
        checkpoint_parameter <- FALSE
      }
      ##
      exECS0037 <- which(PARAM_exECS[, 1] == 'exECS0037')
      if (length(exECS0037) == 0) {
        message("ERROR!!! Problem with exECS0037! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
        checkpoint_parameter <- FALSE
      } else {
        PubChem_library_path <- gsub("\\", "/", PARAM_exECS[exECS0037, 2], fixed = TRUE)
        PARAM_exECS[exECS0037, 2] <- PubChem_library_path
        if (!file.exists(PubChem_library_path)) {
          message("ERROR!!! Problem with exECS0037! PubChem library data is not available! You should use the 'molecular_formula_library_generator' module to produce the molecular formula library!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    if (checkpoint_parameter) {
      message("The spreadsheet is consistent with the IDSL.UFAx workflow!")
      ##
      Elements <- c("As", "Br", "Cl", "Na", "Si", "B", "C", "F", "H", "I", "K", "N", "O", "P", "S") # DO NOT change this order!
      x_c_el <- 7 # = which(Elements == "C") # index number of Carbon
      LElements <- 15  # = length(Elements)
      elementSorterList <- element_sorter(ElementList = Elements, alphabeticalOrder = FALSE)
      Elements_mass_abundance <- elementSorterList[[2]]
      valence_vec <- elementSorterList[[3]]
      Elements_mass <- do.call('c', lapply(Elements, function(el) {
        Elements_mass_abundance[[el]][[1]][1]
      }))
      Elements_mass[6] <- 11.009305167 # To use the most abundant mass of Boron
      ##
      ##########################################################################
      ##
      get_formula_mz <- function(query_mz, qvec, carbon_masses, atom_vec, massAccuracy) {
        ##
        set_1 <- comboGeneral(qvec, 1, FALSE, constraintFun = "sum",
                              comparisonFun = c(">=", "<="),
                              limitConstraints = c((query_mz - massAccuracy), (query_mz + massAccuracy)),
                              keepResults = TRUE)
        
        set_2 <- comboGeneral(qvec, 2, FALSE, constraintFun = "sum",
                              comparisonFun = c(">=", "<="),
                              limitConstraints = c((query_mz - massAccuracy), (query_mz + massAccuracy)),
                              keepResults = TRUE)
        set_2_a <- set_2
        set_2_a[!set_2_a %in% carbon_masses] <- 0
        set_2_a[set_2_a %in% carbon_masses] <- 1
        set_2 <- set_2[which(rowSums(set_2_a) > 0),]
        
        set_3 <- comboGeneral(qvec, 3, FALSE, constraintFun = "sum",
                              comparisonFun = c(">=", "<="),
                              limitConstraints = c((query_mz - massAccuracy), (query_mz + massAccuracy)),
                              keepResults = TRUE)
        set_3_a <- set_3
        set_3_a[!set_3_a %in% carbon_masses] <- 0
        set_3_a[set_3_a %in% carbon_masses] <- 1
        set_3 <- set_3[which(rowSums(set_3_a) > 0),]
        
        set_4 <- comboGeneral(qvec, 4, FALSE, constraintFun = "sum",
                              comparisonFun = c(">=", "<="),
                              limitConstraints = c((query_mz - massAccuracy), (query_mz + massAccuracy)),
                              keepResults = TRUE)
        set_4_a <- set_4
        set_4_a[!set_4_a %in% carbon_masses] <- 0
        set_4_a[set_4_a %in% carbon_masses] <- 1
        set_4 <- set_4[which(rowSums(set_4_a) > 0),]
        
        set_5 <- comboGeneral(qvec, 5, FALSE, constraintFun = "sum",
                              comparisonFun = c(">=", "<="),
                              limitConstraints = c((query_mz - massAccuracy), (query_mz + massAccuracy)),
                              keepResults = TRUE)
        set_5_a <- set_5
        set_5_a[!set_5_a %in% carbon_masses] <- 0
        set_5_a[set_5_a %in% carbon_masses] <- 1
        set_5 <- set_5[which(rowSums(set_5_a) > 0),]
        
        set_6 <- comboGeneral(qvec, 6, FALSE, constraintFun = "sum",
                              comparisonFun = c(">=", "<="),
                              limitConstraints = c((query_mz - massAccuracy), (query_mz + massAccuracy)),
                              keepResults = TRUE)
        set_6_a <- set_6
        set_6_a[!set_6_a %in% carbon_masses] <- 0
        set_6_a[set_6_a %in% carbon_masses] <- 1
        set_6 <- set_6[which(rowSums(set_6_a) > 0), ]
        ##
        get_formula_text <- function(sety) {
          if (!is.matrix(sety)) {
            sety <- t(as.matrix(sety))
          }
          if (nrow(sety) > 0) {
            form_mat <- sapply(1:(ncol(sety) - 1), function(z) {
              as.character(atom_vec[as.character(sety[, z])])
            })
            if (!is.matrix(form_mat)) {
              form_mat <-  t(as.matrix(form_mat))
            }
            apply(form_mat, 1, paste0, collapse = "")
          }
        }
        c(get_formula_text(set_1), get_formula_text(set_2), get_formula_text(set_3),
          get_formula_text(set_4), get_formula_text(set_5), get_formula_text(set_6))
      }
      ##
      ##########################################################################
      ##
      exhaustive_chemical_enumeration_call <- function(i_mz) {
        exECS_annontated_table <- NULL
        ##
        query_mz <- as.numeric(peaklist[i_mz, 8])
        R13C_PL <- as.numeric(peaklist[i_mz, 11])
        ##
        c_min <- max(c(1, floor((R13C_PL - maxR13C)/1.0816)))
        c_max <- ceiling((R13C_PL + maxR13C)/1.0816)
        b_max <- min(c(B_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[6])))))
        br_max <- min(c(Br_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[2])))))
        cl_max <- min(c(Cl_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[3])))))
        k_max <- K_MAX
        s_max <- min(c(S_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[15])))))
        si_max <- min(c(Si_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[5])))))
        n_max <- min(c(N_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[12])))))
        as_max <- min(c(As_MAX, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[1])))
        f_max <- min(c(2*c_max + 3*n_max + 6, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[8])))))
        i_max <- min(c(I_MAX, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[10])))
        h_max <- min(c(2*c_max + 3*n_max + 6, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[9])))))
        na_max <- Na_MAX
        o_max <- min(c(O_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[13])))))
        p_max <- min(c(P_MAX, max(c(0, floor((query_mz - c_min*Elements_mass[x_c_el])/Elements_mass[14])))))
        ##
        minfreq <- rep(0, LElements)
        minfreq[x_c_el] <- c_min        
        maxfreq <- c(as_max, br_max, cl_max, na_max, si_max, b_max, c_max, f_max, h_max, i_max, k_max, n_max, o_max, p_max, s_max)
        #
        qvec <- unlist(lapply(1:LElements, function(x) {
          xvec <- round(sapply(c(minfreq[x]:maxfreq[x]), function(yy) {yy*Elements_mass[x]}), digits = 3)
          names(xvec) <- sapply(c(minfreq[x]:maxfreq[x]), function(yy) {paste0(Elements[x], yy)})
          xvec
        }))
        qvec <- qvec[which(qvec != 0)]
        atom_vec <- names(qvec)
        names(atom_vec) <- as.character(qvec)
        carbon_masses <- Elements_mass[x_c_el]*(minfreq[x_c_el]:maxfreq[x_c_el])
        molecular_formula <- get_formula_mz(query_mz, qvec, carbon_masses, atom_vec, massAccuracy)
        ########################################################################
        if (length(molecular_formula) > 0) {
          MoleFormVecMat <- do.call(rbind, lapply(molecular_formula, function(molf) {
            molvec <- formula_vector_generator(molf, Elements, LElements, allowedRedundantElements = FALSE)
            if (molvec[x_c_el] > 0) {
              molvec
            }
          }))
          ######################################################################
          LMoleFormVecMat <- length(MoleFormVecMat)/15
          if (LMoleFormVecMat > 0) {
            ####################################################################
            ## xCondition <- which((h + cl + br + f + i) >= (c/2 - n - 1) & (h + cl + br + f + i) <= (2*c + 3*n + 6))
            if (LMoleFormVecMat == 1) {
              sumClBrFI <- sum(MoleFormVecMat[c(2, 3, 8, 10)])
              sumHClBrFI <- sumClBrFI + MoleFormVecMat[9]
            } else {
              sumClBrFI <- rowSums(MoleFormVecMat[, c(2, 3, 8, 10)])
              sumHClBrFI <- sumClBrFI + MoleFormVecMat[, 9]
            }
            xCondition <- which(sumHClBrFI >= (1/2*MoleFormVecMat[, x_c_el] - MoleFormVecMat[, 12] - 1) &
                                  sumHClBrFI <= (2*MoleFormVecMat[, x_c_el] + 3*MoleFormVecMat[, 12] + 6))
            sumHClBrFI <- NULL
            ##
            if (length(xCondition) > 0) {
              if (maxHalogenCounts > 0) {
                xC <- which(sumClBrFI[xCondition] <= maxHalogenCounts)
                if (length(xC) > 0) {
                  xCondition <- xCondition[xC]
                } else {
                  xCondition <- NULL
                }
              }
              sumClBrFI <- NULL
              ##
              if (length(xCondition) > 0) {
                if (NaKruleCheck) {
                  xC <- which((MoleFormVecMat[xCondition, 4] + MoleFormVecMat[xCondition, 11]) <= 1)
                  if (length(xC) > 0) {
                    xCondition <- xCondition[xC]
                  } else {
                    xCondition <- NULL
                  }
                }
                ##
                LMoleFormVecMat <- length(xCondition)
                if (LMoleFormVecMat > 0) {
                  if (LMoleFormVecMat == 1) {
                    MoleFormVecMat <- matrix(MoleFormVecMat[xCondition, ], nrow = 1)
                  } else {
                    MoleFormVecMat <- MoleFormVecMat[xCondition, ]
                  }
                  ##
                  seniorRuleMatCheck <- do.call('c', lapply(1:LMoleFormVecMat, function(el) {# Extended SENOIR rule
                    extendedSENIORrule(MoleFormVecMat[el, ], valence_vec, ionization_correction = ipw_n)
                  }))
                  LMoleFormVecMat <- length(which(seniorRuleMatCheck))
                  if (LMoleFormVecMat > 0) {
                    if (LMoleFormVecMat == 1) {
                      MoleFormVecMat <- matrix(MoleFormVecMat[seniorRuleMatCheck, ], nrow = 1)
                    } else {
                      MoleFormVecMat <- MoleFormVecMat[seniorRuleMatCheck, ]
                      MoleFormVecMat <- unique(as.matrix(MoleFormVecMat)) # To remove redundant rows
                    }
                    ##
                    molecularFormulaDatabase <- list(Elements, MoleFormVecMat, rep(1, LMoleFormVecMat))
                    ##
                    MoleFormVecMat <- NULL
                    ##
                    IPDB_exECS <- molecularFormula2IPdb(molecularFormulaDatabase, retentionTime = NULL, peak_spacing, intensity_cutoff_str, IonPathways = "[M]",
                    					number_processing_threads = 1, UFA_IP_memeory_variables, allowedMustRunCalculation = FALSE, allowedVerbose = FALSE)
                    ##
                    exECS_annontated_table <- molecular_formula_annotator(IPDB_exECS, spectraList, peaklist, selectedIPApeaks = i_mz, massAccuracy, maxNEME, minPCS,
                                                                          minNDCS, minRCS, scoreCoefficients, RTtolerance = NA, correctedRTpeaklist = NULL,
                                                                          exportSpectraParameters, number_processing_threads = 1)
                  }
                }
              }
            }
          }
        }
        save(exECS_annontated_table, file = paste0(output_path_log, i_mz, '_exECS_annotation.Rdata'))
        ##
        return(exECS_annontated_table)
      }
      ##########################################################################
      outputer <- IDSL.IPA::IPA_MSdeconvoluter(inputHRMSfolderPath, MSfileName)
      spectraList <- outputer[["spectraList"]]
      msPolarity <- outputer[["msPolarity"]]
      outputer <- NULL
      ## Creating the log folder
      output_path_log <- paste0(output_path, "/log_exECS_annotation_", MSfileName, "/")
      if (!dir.exists(output_path_log)) {
        tryCatch(dir.create(output_path_log, recursive = TRUE), warning = function(w){stop("Can't create the log folder!!!")})
      }
      ##
      if (exportSpectraCheck) {
        exportSpectraParameters <- c(maxAllowedNumberHits, MSfileName, msPolarity, outputProfileSpectra)
      } else {
        exportSpectraParameters <- NULL
      }
      ##
      message("Initiated the exhaustive chemical enumeration analysis!!!")
      ##
      ##########################################################################
      ##
      if (number_processing_threads == 1) {
        ##
        exhaustive_chemical_enumeration_annotated_table <- do.call(rbind, lapply(selectedIDpeaklist, function(i_mz) {
          exhaustive_chemical_enumeration_call(i_mz)
        }))
        ##
      } else {
        ##
        ########################################################################
        ##
        osType <- Sys.info()[['sysname']]
        if (osType == "Windows") {
          clust <- makeCluster(number_processing_threads)
          registerDoParallel(clust)
          ##
          exhaustive_chemical_enumeration_annotated_table <- foreach(i_mz = selectedIDpeaklist, .combine = 'rbind', .verbose = FALSE) %dopar% {
            exhaustive_chemical_enumeration_call(i_mz)
          }
          ##
          stopCluster(clust)
          ##
          ######################################################################
          ##
        } else if (osType == "Linux") {
          ##
          exhaustive_chemical_enumeration_annotated_table <- do.call(rbind, mclapply(selectedIDpeaklist, function(i_mz) {
            exhaustive_chemical_enumeration_call(i_mz)
          }, mc.cores = number_processing_threads))
          closeAllConnections()
        }
      }
      ##
      ##########################################################################
      ##
      rownames(exhaustive_chemical_enumeration_annotated_table) <- NULL
      if (tolower(exECS0035) == "yes") {
        message("Initiated searching in the library of molecular formula!!!")
        MFlibrary <- IDSL.IPA::loadRdata(PubChem_library_path)
        exhaustive_chemical_enumeration_annotated_table <- molecular_formula_library_search(exhaustive_chemical_enumeration_annotated_table, MFlibrary, IonPathways, number_processing_threads)
        message("Completed searching in the library of molecular formula!!!")
      }
      ##
      save(exhaustive_chemical_enumeration_annotated_table, file = paste0(output_path, "/exhaustive_chemical_enumeration_annotated_table_", MSfileName, ".Rdata"))
      write.csv(exhaustive_chemical_enumeration_annotated_table, file = paste0(output_path, "/exhaustive_chemical_enumeration_annotated_table_", MSfileName, ".csv"), row.names = TRUE)
      message("Completed the exhaustive chemical enumeration analysis!!!")
      required_time <- Sys.time() - initiation_time
      print(required_time)
    }
  }
  return(exhaustive_chemical_enumeration_annotated_table)
}