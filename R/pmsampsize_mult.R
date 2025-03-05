pmsampsize_mult_general <- function(type,
                                    nagrsquared = NA,
                                    csrsquared = NA,
                                    rsquared = NA,
                                    parameters,
                                    shrinkage = 0.9,
                                    prevalence = NA,
                                    cstatistic = NA,
                                    seed = 123456,
                                    rate = NA,
                                    timepoint = NA,
                                    meanfup = NA,
                                    intercept = NA,
                                    sd = NA,
                                    mmoe=1.1,
                                    ### New parameters
                                    K = NA,
                                    mult_n_events = NA,
                                    mult_nagrsquared_overall = NA,
                                    mult_rsquared_overall = NA) {

  ### NB: Unless specified as "_overall", metrics (nagrsquared, csrsquared, cstatistic) are referring to pairwise

  ### Run error check
  pmsampsize_errorcheck_mult(type = type,
                             nagrsquared = nagrsquared,
                             csrsquared = csrsquared,
                             rsquared = rsquared,
                             parameters = parameters,
                             shrinkage = shrinkage,
                             cstatistic = cstatistic,
                             prevalence = prevalence,
                             rate = rate,
                             timepoint = timepoint,
                             meanfup = meanfup,
                             intercept = intercept,
                             sd = sd,
                             mmoe=1.1,
                             ### New parameters
                             K = K,
                             mult_n_events = mult_n_events,
                             mult_nagrsquared_overall = mult_nagrsquared_overall,
                             mult_rsquared_overall = mult_rsquared_overall)

  if (type == "m"){

    ### Define number of pairs
    n_pairs <- as.numeric(ncol(combn(K,2)))

    ### Calculate p_k (proportion of individuals in outcome category k)
    p_k <- mult_n_events/sum(mult_n_events)

    ### Calculate p_k_r, proportion of individuals in outcome category k and r combined
    p_k_r <- lapply(1:n_pairs, function(x){
      ### Get locations
      k <- combn(K,2)[1,x]
      j <- combn(K,2)[2,x]

      ### Assign p_k_r
      return((mult_n_events[j] + mult_n_events[k])/sum(mult_n_events))

    })
    p_k_r <- unlist(p_k_r)

    ### Calculate pairwise outcome proportions (phi), of category k relative to category r
    phi <- lapply(1:n_pairs, function(x){
      ### Get locations
      k <- combn(K,2)[1,x]
      j <- combn(K,2)[2,x]

      ### Assign p_k_r
      return((mult_n_events[j])/(mult_n_events[j] + mult_n_events[k]))

    })
    phi <- unlist(phi)

    ### Calculate csrsquared for each pairwise model
    if (any(is.na(csrsquared))){
      if (any(is.na(cstatistic))){
        if (any(is.na(nagrsquared))){
          stop("One of crsquared, nagsquared or cstatistic must be specified")
        } else {

          ### Estimate max_rsquared
          max_rsquared <- lapply(1:n_pairs, function(x) {
            lnLnull_n <- phi[x]*log(phi[x]) + (1-phi[x])*log(1-phi[x])
            max_rsquared <- (1 - exp(2*lnLnull_n))
            return(max_rsquared)
          })

          ### Get max_rsquared
          max_rsquared <- unlist(max_rsquared)

          ### Get csrsquared
          csrsquared <- nagrsquared*max_rsquared

        }
      } else {

        ### Estimate csrsquared from pairwise cstatistics
        approx_rsq <- lapply(1:n_pairs, function(x) {
          approx_rsq_interior <- pmsampsize:::cstat2rsq(cstatistic=cstatistic[x], prevalence=phi[x], seed=seed)
          return(approx_rsq_interior$R2.coxsnell)
        })

        csrsquared <- unlist(approx_rsq)
      }
    } else {

    }

    ### If nagrsquared not specified, calculate pairwise max_rsquared and nagrsquared
    if (any(is.na(nagrsquared))){
      ### Estimate max_rsquared
      max_rsquared <- lapply(1:n_pairs, function(x) {
        lnLnull_n <- phi[x]*log(phi[x]) + (1-phi[x])*log(1-phi[x])
        max_rsquared <- (1 - exp(2*lnLnull_n))
        return(max_rsquared)
      })

      ### Get max_rsquared
      max_rsquared <- unlist(max_rsquared)
      print(max_rsquared)

      ### Get nagelkerke rsquared
      nagrsquared <- csrsquared/max_rsquared
    }

    ### Calculate max(R2_CS) for the overall model
    max_mult_rsquared_overall <- 1 - prod(p_k^p_k)^2

    ### Calculate R2_CS_adj for the overall model
    if (is.na(mult_rsquared_overall)){

      ### Calculate an estimte of R2_CS_app, based off R2_NAGEL = 0.15
      mult_rsquared_overall <- mult_nagrsquared_overall*max_mult_rsquared_overall
    }

    ###########################
    ### Create output object ##
    ###########################
    out <- pmsampsize_mult(csrsquared=csrsquared,
                           parameters=parameters,
                           shrinkage=shrinkage,
                           cstatistic=cstatistic,
                           nagrsquared=nagrsquared,
                           max_rsquared = max_rsquared,
                           K = K,
                           n_pairs = n_pairs,
                           p_k = p_k,
                           p_k_r = p_k_r,
                           phi = phi,
                           mult_rsquared_overall = mult_rsquared_overall,
                           max_mult_rsquared_overall = max_mult_rsquared_overall)

  }

  return(out)

}


pmsampsize_mult <- function(csrsquared,parameters,shrinkage,cstatistic,nagrsquared, max_rsquared,
                            K, n_pairs, p_k, p_k_r, phi, mult_rsquared_overall, max_mult_rsquared_overall) {

  ####################
  ### Criterion i) ###
  ####################

  ### Calculate m_k_r using pmsampsize_bin
  m_k_r <- lapply(1:n_pairs, function(x) {
    pmsampsize:::pmsampsize_bin(rsquared = csrsquared[x],
                                parameters = parameters,
                                prevalence = phi[x],
                                shrinkage = 0.9,
                                cstatistic = cstatistic[x])$results_table["Criteria 1", "Samp_size"]
  })
  m_k_r <- unlist(m_k_r)

  ### Calculate n_k_r for criterion (i) for each submodel
  n_k_r <- m_k_r/p_k_r

  ### Take the ceiling of the maximum of these as the sample size for criteiron (i)
  N_C1 <- ceiling(max(n_k_r))

  #####################
  ### Criterion ii) ###
  #####################

  N_C2 <- (K-1)*parameters/((mult_rsquared_overall/(mult_rsquared_overall + 0.05*max_mult_rsquared_overall) - 1)*log(1 - mult_rsquared_overall - 0.05*max_mult_rsquared_overall))
  N_C2 <- ceiling(N_C2)

  ######################
  ### Criterion iii) ###
  ######################

  N_C3_vec <- lapply(1:K, function(x) {qchisq(1-0.05/5, 1)*p_k[x]*(1-p_k[x])/0.05^2})
  N_C3_vec <- unlist(N_C3_vec)

  N_C3 <- ceiling(max(N_C3_vec))

  ###########################
  ### create output table ###
  ###########################

  # minimum n
  nfinal <- max(N_C1,N_C2,N_C3)

  #######################################################################
  ### Calculate expected shrinkage for each pair at given sample size ###
  #######################################################################
  shrinkage_targeted <- lapply(1:n_pairs, function(x) {1 + parameters/(nfinal*p_k_r[x]*log(1 - csrsquared[x]/shrinkage))})
  shrinkage_targeted <- unlist(shrinkage_targeted)

  ###########################
  ### Create output object ##
  ###########################

  #########
  ### Make a table for criteria 1
  #########

  # create output table
  res_criteria1 <- matrix(NA,n_pairs, 4)
  colnames(res_criteria1) <- c("CS_Rsq", "Max_Rsq","Nag_Rsq", "targ_shrinkage")
  rownames(res_criteria1) <- paste("Pair", unlist(lapply(1:n_pairs, function(x){paste(combn(K,2)[,x], collapse = ",")})), sep = " ")
  res_criteria1[,1] <- round(csrsquared, 3)
  res_criteria1[,2] <- round(max_rsquared, 3)
  res_criteria1[,3] <- round(nagrsquared, 3)
  res_criteria1[,4] <- round(shrinkage_targeted, 3)

  #########
  ### Make a table for criteria 2
  #########

  #########
  ### Make a table for criteria 3
  #########


  #########
  ### Make overall table
  #########
  res_final <- c(N_C1,N_C2,N_C3,nfinal)
  names(res_final) <- c("Criteria 1","Criteria 2","Criteria 3","Final")

  output_object <- list(results_criteria1 = res_criteria1,
                        results_all = res_final,
                        sample_size = nfinal,
                        parameters = parameters,
                        events_expected = nfinal*p_k,
                        type = "multinomial")

  return(output_object)

}

pmsampsize_errorcheck_mult <- function(type,
                                       csrsquared,
                                       nagrsquared,
                                       rsquared,
                                       parameters,
                                       shrinkage,
                                       cstatistic,
                                       prevalence,
                                       rate,
                                       timepoint,
                                       meanfup,
                                       intercept,
                                       sd,
                                       mmoe,
                                       ### New parameters
                                       K,
                                       mult_n_events,
                                       mult_nagrsquared_overall,
                                       mult_rsquared_overall){

  if (type == "m"){

    ### Write NA test that works when alternative could be vectors
    my_na_test <- function(input){isTRUE(length(input) == 1 & is.na(input))}

    # Ensure K and mult_n_events specified
    if (is.na(K)){
      stop("K must be specified")
    }
    if (my_na_test(mult_n_events)){
      stop("mult_n_events must be specified")
    }

    # Make sure prevalence not specified
    if (!is.na(prevalence)){
      stop("Prevalence should not be specified. Prevalence is estimated from number of events expected in each outcome category, specified through mult_n_events")
    }

    ### Define the number of events in each category
    if (length(mult_n_events) != K){
      stop("The length of mult_n_events should be equal to the number of outcome categories.")
    }

    ### Only specify one of csrsquared, cstatistic, nagrsquared
    if (sum(as.numeric(!my_na_test(cstatistic)), as.numeric(!my_na_test(csrsquared)), as.numeric(!my_na_test(nagrsquared))) > 1){
      stop("Specify only one of csrsquared, nagrsquared, or C-statistic as a numeric vector with no NA values")
    }
    ### Specify one of csrsquared, cstatistic, nagrsquared
    if (sum(as.numeric(!my_na_test(cstatistic)), as.numeric(!my_na_test(csrsquared)), as.numeric(!my_na_test(nagrsquared))) == 0){
      stop("Specify one of csrsquared, nagrsquared, or C-statistic as a numeric vector with no NA values")
    }


    ### check not all specified
    if (any(is.na(csrsquared))){
      if (any(is.na(cstatistic))){
        if (any(is.na(nagrsquared))){
          stop("Specify one of csrsquared, nagrsquared, or C-statistic as a numeric vector with no NA values")
        } else {
          if (length(nagrsquared) != ncol(combn(K,2))){
            stop("The length of nagrsquared should be equal to the number of pairs of outcome categories. This is equal to K choose 2")
          }

        }
      } else {

        if (length(cstatistic) != ncol(combn(K,2))){
          stop("The length of cstatistic should be equal to the number of pairs of outcome categories. This is equal to K choose 2")
        }

        if (!any(is.na(nagrsquared))){
          stop("Specify only one of csrsquared, nagrsquared, or C-statistic")
        }
      }
    } else {
      ## CAN THIS BE ADDED IN TO ERROR CHECK???
      if (!any(is.na(nagrsquared))){
        stop("Specify only one of csrsquared, nagrsquared, or C-statistic")
      }
      if (!any(is.na(cstatistic))){
        stop("Specify only one of csrsquared, nagrsquared, or C-statistic")
      }
      if (length(csrsquared) != ncol(combn(K,2))){
        stop("The length of csrsquared should be equal to the number of pairs of outcome categories. This is equal to K choose 2")
      }
    }

    ### cstatistic
    if (all(!is.na(cstatistic))){
      if (is.numeric(cstatistic) == FALSE){
        stop("cstatistic values must be numeric")
      }
      if (any(cstatistic < 0 | cstatistic > 1)){
        stop("cstatistic values must be between 0 and 1")
      }
    }

    ### csrsquared
    if (all(!is.na(csrsquared))){
      if (is.numeric(csrsquared) == FALSE){
        stop("csrsquared values must be numeric")
      }
      if (any(csrsquared < 0 | csrsquared > 1)){
        stop("csrsquared values must be between 0 and 1")
      }
    }

    ### nagrsquared
    if (all(!is.na(nagrsquared))){
      if (is.numeric(nagrsquared) == FALSE){
        stop("nagrsquared values must be numeric")
      }
      if (any(nagrsquared < 0 | nagrsquared > 1)){
        stop("nagrsquared values must be between 0 and 1")
      }
    }

    ### Overall model metrics
    if (my_na_test(mult_rsquared_overall)){
      if (my_na_test(mult_nagrsquared_overall)){
        stop("Specify one of mult_rsquared_overall or mult_nagrsquared_overall")
      }
      if (length(mult_nagrsquared_overall) > 1){
        stop("mult_nagrsquared_overall should be a numeric of length 1")
      }
      if (!is.numeric(mult_nagrsquared_overall)){
        stop("mult_rsquared_overall should be a numeric of length 1")
      }
    } else {

      if (!my_na_test(mult_nagrsquared_overall)){
        stop("Specify only one of mult_nagrsquared_overall or mult_rsquared_overall")
      }

      if (length(mult_rsquared_overall) > 1){
        stop("mult_rsquared_overall should be a numeric of length 1")
      }
      if (!is.numeric(mult_rsquared_overall)){
        stop("mult_rsquared_overall should be a numeric of length 1")
      }
    }
    ### TESTS OVER

  }
}


