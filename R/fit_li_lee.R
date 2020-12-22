#' Fit a multipopulation Li-Lee mortality model
#'
#' @description This function fits a multipopulation Li-Lee model. This model consists of two times a
#'   Lee & Carter model: one for the common trend and one for the country-specific
#'   deviation from this trend. The model has the following form
#'   \describe{
#'     \item{}{\eqn{log(mtx) = log(mtx^{(T)}) + log(mtx^{(c)})}}
#'     \item{}{\eqn{log(mtx^{(T)}) = Ax + Bx * Kt}}
#'     \item{}{\eqn{log(mtx^{(c)}) = ax + bx * kt},}
#'     }
#'     where \eqn{log(mxt)} refers to the force of mortality for an x-year old at time t.
#'
#' @param xv The vector of ages.
#' @param yv The vector of years.
#' @param yv_spec The vector of years for the country of interest.
#' @param country_spec The country of interest.
#' @param data The mortality data (dtx, wtx, wa) for each country.
#' @param method The method to fit the mortality model (currently only \code{"NR"} working).
#' @param Ax Logical. Should a common age specific factor be included in the model?
#' @param exclude_coh Logical. If TRUE, the weights are set to zero for cohorts with
#'   fewer than 5 observations.
#'
#' @details The \code{data} must be (in the same format as) the output of the function
#'   \code{\link[MultiMoMo]{get_mortality_data}}.
#'
#'   The parameters are fitted using a two step Newton-Raphson approach.
#'   In the first step, we fit the parameters for the commond trend
#'   (\eqn{Ax}, \eqn{Bx} and \eqn{Kt}) (standard Lee-Carter fit).
#'   In a second step, we keep \eqn{Ax}, \eqn{Bx} and \eqn{Kt} fixed and assume
#'   that \deqn{dtx ~ POI(etx*mtx^{(T)}*mtx^{(c)}),} where \eqn{mtx^{(T)}}
#'   denote the fitted forces of mortality for the common trend, namely
#'   \deqn{mtx^{(T)} = exp(Ax + Bx*Kt),}
#'   and where the \eqn{mtx^{(c)}} denote the fitted forces of mortality
#'   for the country-specific deviation from the common trend, namely
#'   \deqn{mtx^{(c)} = exp(ax + bx* kt).}
#'   So in the second step, we again fit a standard Lee-Carter model for the deaths
#'   of your country of interest \eqn{dtx} and with the adapted exposures \eqn{etx*mtx^{(T)}}.
#'
#'
#'   The BIC is defined as:
#'   \deqn{BIC = -2*ll + log(sum(wa)).}
#'
#'   The log-likelihood \eqn{ll} equals:
#'   \deqn{ll = sum((dtx*log(etx*mtx) - etx*mtx - lgamma(dtx + 1))*wa).}
#'   Here, \eqn{mtx} denotes the fitted forces of mortality for your country of interest,
#'   \eqn{etx} its exposure numbers and \eqn{dtx} the death counts.
#'
#'   The Poisson Pearson residuals equal:
#'   \deqn{(dtx - etx*mtx)/sqrt(etx*mtx).}
#'
#'
#' @return A list containing the following objects
#' \itemize{
#'   \item Common parameters: $Ax, $Bx and $Kt
#'   \item Country specific parameters: $ax, $bx and $kt
#'   \item The force of mortality estimates: $mtx
#'   \item The log-likelihood of the model: $ll
#'   \item The number of parameters in the model: $npar
#'   \item The BIC of the model: $BIC
#'   \item The Poisson Pearson residuals of the mortality model: $residuals
#' }
#'
#' @examples
#' lst   <- MultiMoMo::european_mortality_data
#' dat_M <- lst$M
#' dat_F <- lst$F
#' xv    <- 0:90
#' yv = yv_spec <- 1970:2018
#' Countries   <- names(dat_M$UNI)
#' country_spec <- "BE"
#' fit_M <- fit_li_lee(xv, yv, yv_spec, country_spec, dat_M, "NR", TRUE, FALSE)
#' fit_F <- fit_li_lee(xv, yv, yv_spec, country_spec, dat_F, "NR", TRUE, FALSE)
#'
#' @export


fit_li_lee <- function(xv, yv, yv_spec, country_spec, data, method, Ax, exclude_coh){

  if(method != "NR")
    stop("Only the Newton-Raphson method (NR) works for now.")

  # First, a Lee-Carter model for the common evolution of mortality is constructed
  dtxALL <- data$ALL$dtx
  etxALL <- data$ALL$etx
  waALL  <- data$ALL$wa

  LC.total <- fit701(xv, yv, etxALL, dtxALL, waALL, xv*0, "ALL", exclude_coh)

  A.x <- LC.total$beta1
  B.x <- LC.total$beta2	   	# note: sum B.x^2 = 1
  K.t <- LC.total$kappa2 		# note: sum K.t = 0

  # Project K.t for difference in data years between all countries and the country of interest
  l = tail(yv,1)
  diff = tail(yv_spec,1) - l
  if(diff>0){
    K.t.proj = K.t[length(yv)] + (1:diff)*(K.t[length(yv)]-K.t[1])/(yv[length(yv)]-yv[1])
    K.t = c(K.t,K.t.proj)
  }

  # set A.x = 0 if Ax = FALSE
  if(Ax == FALSE)
    A.x <- A.x*0

  # Second a Lee-Carter model for the country-specific deviation from the overall trend
  alpha.x	 <- c()	  # initialize list alphas
  beta.x	 <- c()	  # initialize list betas
  kappa.t	 <- c()	  # initialize list kappas
  m.txi	   <- c() 	# initialize list death rates

  c <- which(names(data$UNI) == country_spec)
  data$UNI[[c]]$dtx[which(data$UNI[[c]]$dtx == 0)] <- 1e-12

  # alpha, beta, kappa for BE
  LC <-  fit701(xv, yv_spec, (data$UNI[[c]]$etx[as.character(yv_spec),]) * exp(rep(A.x, each = length(yv_spec), times = 1)
                                                      + tail(K.t, length(yv_spec)) %o% B.x),
                data$UNI[[c]]$dtx[as.character(yv_spec),], data$UNI[[c]]$wa[as.character(yv_spec),],
                xv*0, "ALL", exclude_coh) # adjust exposure!

  alpha.x	  <- LC$beta1
  beta.x 	  <- LC$beta2		   	# note: sum = 1
  kappa.t   <- LC$kappa2	 		# note: sum = 0
  ll        <- LC$ll
  npar      <- LC$npar + LC.total$npar
  BIC       <- -2*LC$ll + log(sum(LC$wa))*npar
  residuals <- LC$epsilon

  # death rates country_spec
  m.tx = t(exp(rep(A.x, each = 1, times = length(yv_spec)) +
                 rep(alpha.x, each = 1, times = length(yv_spec)) +
                 B.x%o%tail(K.t, length(yv_spec)) + beta.x%o%kappa.t))
  dimnames(m.tx) <- list(yv_spec, xv)

  results <- list(A.x = A.x, B.x = B.x, K.t = K.t, a.x = alpha.x, b.x = beta.x, k.t = kappa.t,
              m.tx = m.tx, ll = ll, npar = npar, BIC = BIC, residuals = residuals)
  results
}


# LifeMetrics Software
# License Agreement
#
# This software license agreement ("Agreement") is an agreement between JPMorgan Chase Bank, N.A.,
# (collectively with its affiliates, "Morgan") and you (either an individual or an entity on whose behalf the individual is acting)
# ("you" or "Licensee") regarding your use of the LifeMetrics Software by JPMorgan and is effective on the date that you
# either compile the Software Source Code or otherwise use the Products ("Agreement Effective Date").
#
# BEFORE YOU EITHER COMPILE THE SOFTWARE SOURCE CODE OR OTHERWISE USE THE PRODUCTS, CAREFULLY READ THE TERMS AND CONDITIONS
# OF THIS AGREEMENT BELOW. BY COMPILING THE SOFTWARE SOURCE CODE OR USING THE PRODUCTS IN ANY MANNER, YOU ARE AGREEING TO BE
# BOUND BY AND ARE BECOMING A PARTY TO THIS AGREEMENT. IF YOU DO NOT AGREE TO ALL OF THE TERMS OF THIS AGREEMENT,
# IMMEDIATELY DELETE ALL COPIES OF THE PRODUCTS FROM YOUR COMPUTER SYSTEMS.
#
#
# 1. Grant of License.
#
# (a) Subject to the terms and conditions of this Agreement, Morgan hereby grants to Licensee, a non-transferable,
# non-assignable, personal, revocable, non-exclusive license: (i) to compile the source code of the LifeMetrics software
# including all components and updates thereto, and to the proprietary data contained therein, (the "Software Source Code")
# using a third party product, the results of which is the Software in object code form, ("Software Object Code"),
# ("Software Source Code" and the "Software Object Code" are hereinafter collectively referred to as the "Software"),
# and (ii) to use the Software to evaluate mortality projections, and (iii) to use any Software user manuals provided by Morgan,
# ("Manuals").  The "Software" and the "Manuals" are hereinafter referred to as "Products".
#
# (b) No rights to use the Products are granted hereunder other than those expressly granted in Section 1(a).
# Except as expressly authorized in this Agreement, Licensee shall not rent, lease, sublicense, distribute, transfer,
# copy, reproduce, display, modify or timeshare the Products or any portion thereof, or use the Products as a
# component of or a base for products or services prepared for commercial or non commercial sale, sublicense, lease,
# access, hosting, application service providing distribution outside the Licensee's organization, or prepare any derivative work
#  based on the Products. Licensee shall not allow any third party or unlicensed user or computer system to access or use the Products.
# Licensee shall not reverse engineer, modify in any way, or create derivative works from the Products, or any portion thereof.
# Licensee agrees that all improvements and modifications to the Products or any part thereof (whether developed by Morgan,
# Licensee or any third party acting on behalf of them at any time shall be and remain the sole and exclusive property of Morgan.
#
# (c) A hyperlink to a software product provided by a third party that may be used to compile the Software Source Code
# may be provided by Morgan ("Third Party Product").  Morgan does not endorse or assume any responsibility or liability
# with regards to the Third Party Product.
#
# 2. Term; Termination.
#
# This Agreement commences on the Agreement Effective Date and expires twelve (12) months thereafter ("Term").
# Morgan may terminate this Agreement at any time for any reason upon thirty (30) days notice, provided,
# however if Licensee breaches any material term of this Agreement, Morgan may terminate this Agreement immediately without any notice.
#
# 3. License fee.
#
# There is currently no license fee for the use of the Products hereunder.
#
# 4. Ownership.
#
# Morgan retains all right, tile and interest to the Products (including any corrections, updates, adaptations,
# enhancements or copies).  The Products are the confidential information of Morgan.  Licensee agrees to take all reasonable
# steps to protect the Products from unauthorized copying or use.
#
# 5. Representations, Warranties, Disclaimers, Indemnification.
#
# (a) Morgan and Licensee each represent and warrant as to itself that (i) it is duly organized, validly existing,
# and in good standing under the laws of the jurisdiction of its organization or incorporation, (ii) it has the power to execute
# the Agreement and to perform its obligations under the Agreement, and (iii) such execution and performance do not violate
# or conflict with any law applicable to it,, any order or judgment of any court or other agency of government applicable to it
# or any of its assets, or any contractual restriction binding on or affecting it or any of its assets.
#
# (b) EXCEPT AS SET FORTH IN THIS SECTION 6, MORGAN MAKES NO WARRANTY, REPRESENTATION, CONDITION, OR AGREEMENT
# WITH RESPECT TO THE PRODUCT OR THE THIRD PARTY PRODUCT.  MORGAN EXPRESSLY DISCLAIMS AND EXCLUDES TO THE FULLEST EXTENT PERMITTED
# BY APPLICABLE LAW ANY AND ALL IMPLIED WARRANTIES, INCLUDING ANY IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE
# FOR THE PRODUCT AND THE THIRD PARTY PRODUCT.  MORGAN SHALL HAVE NO LIABILITY, CONTINGENT OR OTHERWISE, TO LICENSEE OR TO THIRD PARTIES,
# FOR THE ACCURACY, COMPLETENESS OR CURRENCY OF THE PRODUCT OR THE THIRD PARTY PRODUCT OR FOR DELAYS OR OMISSIONS THEREIN, OR FOR INTERRUPTIONS
# IN THE DELIVERY OF THE PRODUCT OR THE THIRD PARTY PRODUCT.  IN NO EVENT WILL MORGAN BE LIABLE FOR ANY DIRECT, ORDINARY,
# SPECIAL, INDIRECT, INCIDENTAL AND/OR CONSEQUENTIAL DAMAGES WHICH MAY BE INCURRED OR EXPERIENCED BY LICENSEE AND/OR ANY
# THIRD PARTIES INCLUDING BUT NOT LIMITED TO ON ACCOUNT OF LICENSEE ENTERING INTO AND/OR LICENSEE AND/OR ANY THIRD PARTY
# RELYING ON THIS AGREEMENT, EVEN IF MORGAN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
#
# (c) Licensee shall indemnify, protect, and hold harmless Morgan, its affiliates, employees and agents,
# from and against any and all losses, liabilities, judgments, suits, actions, proceedings, claims, damages, costs
# (including reasonable attorney's fees) resulting from or arising out of the use by Licensee of the Products and/or the Third Party Product.
#
# 6. Injunctive Relief.
#
# Licensee acknowledges that a breach of any provision of this Agreement will cause Morgan party irreparable
# injury and damage and therefore Licensee may be enjoined through injunctive proceedings in addition to any other rights
# and remedies which may be available to Morgan at law or in equity without the posting of a bond.
#
# 7. General.
#
# (a) All notices and other communications under this Agreement shall be in writing and deemed given upon receipt.
# Notwithstanding the foregoing, Morgan may provide Licensee with notice by posting such Notice it Morgan's Web site and
# such notice shall be deemed valid in accordance with this Agreement.  Notice to Morgan shall be made to:
# JPMorgan Chase Bank, N.A., Technology, Sourcing, and IP Practice Group, 25th Floor, 1 Chase Manhattan Plaza, New York, NY 10081, Fax:  212-383-0800.
#
# (b) The following sections shall survive termination of this Agreement: 4, 5, 6, and 7.
#
# (c) The Agreement constitutes the entire agreement and understanding of the parties with
# respect to its subject matter and supersedes all oral communications and prior writings with respect thereto.
# Morgan may revise this Agreement by posting such revisions to its Web site and your continued use of the Products
# shall constitute your acceptance of such new terms.
#
# (d) Licensee cannot assign this Agreement nor any interest or obligation in or under the Agreement without Morgan's prior written consent.
#
# (e) Morgan will not be liable for any failure, delay or omission to perform any of its obligations
# under the Agreement arising from any cause beyond its reasonable control, including, without limitation, acts of God,
# acts or regulations of government or other authorities, war, fire, strikes or other industrial disputes, power failure,
# failure of telecommunication lines, connection or equipment, or failure or defects in any hardware or software owned or supplied by third parties.
#
# (f) If any provision of the Agreement is held invalid, illegal or unenforceable, the validity,
# legality or enforceability of the remainder of the Agreement will not in any way be affected or impaired thereby.
#
# (g) Licensee hereby agrees not to export or transmit the Products to any country to which such
# transmission is restricted by applicable regulations or statutes.  Licensee agrees to permit Morgan access to
# Licensee's location where the Products is installed and/or used solely to permit Morgan to ensure compliance with the terms of this Agreement.
#
# (h) The Agreement will be governed by the laws of the State of New York without reference to choice of law principles.
# Any action or proceeding relating to this Agreement may be brought and enforced in the courts of the State of New York
# or the United States District Court located in the Borough of Manhattan in New York City. Each of the parties
# irrevocably submits to the jurisdiction of such courts in connection with any such action or proceeding.
# Any process or other legal summons in connection with any such action or proceeding may be served by mailing
# a copy thereof by certified or registered mail, or any substantially similar form of mail, addressed to a party as provided
# for notices hereunder. THE PARTIES WAIVE TRIAL BY JURY IN RESPECT OF ANY PROCEEDINGS RELATED TO THE AGREEMENT.



#' Newton-Raphson implementation to fit a single Lee-Carter mortality model (M1)
#'
#' @description This internal function fits a Lee-Carter mortality model (type M1) with the use
#'   of Newton-Raphson. This mortality model M1 has the following form:
#'   \deqn{log m(t,x) = beta1(x) + beta2(x)*kappa2(t) + Poisson error}
#'
#' @param xv The vector of ages.
#' @param yv The vector of years.
#' @param etx m x n matrix of exposures.
#' @param dtx m x n matrix of death counts.
#' @param wa  m x n matrix of weights (0 or 1).
#' @param int The prespecified vector for \eqn{beta1(x)}. When the zero-vector is supplied,
#'   the Newton_Raphson algorithm calibrates \eqn{beta1(x)}
#' @param constraints Character string. If \code{constraints = "ALL"}, two
#'   identifiability constraints are imposed (see details). If
#'   \code{constraints == "ADJUST"}, only the constraint on beta2(x) is used.
#' @param exclude_coh logical. If TRUE, the weights are set to zero for cohorts with
#'   fewer than 5 observations.
#'
#' @details The two identifiability constraints are:
#'   \deqn{kappa2(t1) = 0,    sum(beta2(x)^2) = 1.}
#'   The BIC is defined as:
#'   \deqn{-2*l1 + log(sum(wa)),}
#'   and where the log-likelihood equals:
#'   \deqn{sum((dtx*log(etx*mhat.tx) - etx*mhat.tx - lgamma(dtx + 1))*wa).}
#'   Here, \eqn{mhat.tx} denotes the fitted forces of mortality, \eqn{etx} the
#'   exposure numbers and \eqn{dtx} the death counts.
#'
#' @return A list containing the fitted parameter estimates of the type M1 Lee-Carter model, as
#'   well as some model specifications.
#'
#'
#'
#' @keywords internal

fit701 <- function(xv, yv, etx, dtx, wa, int, constraints, exclude_coh){

  mtx	<- dtx/etx      			# matrix of death rates
  qtx	<- 1-exp(-mtx)   			# matrix of mortality rates

  n   <- length(xv)		  		# number of ages
  m   <- length(yv)		  		# number of years

  cy  <- (yv[1]-xv[n]):(yv[m]-xv[1])  	# cohort approximate years of birth

  # Initialise parameter vectors for the NR algorithm
  beta1v  <- int
  beta2v  <- (1:n)*0
  beta3v  <- (1:n)*0			          # dummy vector, this will stay at 0
  kappa2v <- (1:m)*0
  gamma3v <- (1:(n+m-1))*0      		# dummy vector, this will stay at 0
  ia      <- array((1:m),c(m,n))  	# matrix of year indexes, i, for the data
  ja      <- t(array((1:n),c(n,m)))	# matrix of age indexes, j, for the data
  ya      <- ia-ja		 	          	# matrix of year of birth indexes for the data
  imj     <- (1-n):(m-1)	      		# the range of values taken by i-j
  lg      <- n+m-1		 		          # number of different values taken by i-j
  ca      <- ya+yv[1]-xv[1]	       	# matrix of years of birth

  if(exclude_coh == TRUE){
    for(k in 1:lg){
      nk <- sum((ca == cy[k])*wa)
      if(nk < 5)
        wa <- wa*(1- (ca == cy[k]))}}

  ww <- cy*0+1	 # this is a vector of 1's and 0's with a 0 if the cohort is completely excluded


  # Stage 0
  # Gives initial estimates for beta1(x), beta2(x) and kappa2(t)
  mx <- mean(xv)
  for(j in 1:n){
    if(all(int == 0))
      beta1v[j] <- sum(log(mtx[,j])*wa[,j])/sum(wa[,j])
    beta2v[j] <- 1/n}

  kappa2v <- (m:1) - (m+1)/2

  # Stage 1: iterate
  l0 <- -1000000
  l1 <- -999999
  iteration <- 0

  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001

  while(abs(l1 - l0) > 1e-12)
  {
    iteration <- iteration + 1

    l0 <- l1

    # Stage 1B optimise over the beta2(x)
    for(j in 1:n){
      # cycle through the range of ages
      dv <- dtx[,j]	# actual deaths
      ev <- etx[,j]	# exposure
      beta2v[j] <- llmaxM2B(beta1v[j], beta2v[j], beta3v[j],
                         kappa2v, gamma3v[(n+1-j):(n+m-j)], dv, ev, wv = wa[,j])
    }

    mhat <- mtx*0
    for(i in 1:m)
      mhat[i,] <- exp(beta1v + beta2v*kappa2v[i] + beta3v*gamma3v[(n+i-1):i])

    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat)
    l1      <- sum((dtx*log(etx*mhat) - etx*mhat - lgamma(dtx + 1))*wa)

    # Stage 1D optimise over the kappa2(t)
    for(i in 1:m){
      # cycle through the range of years
      dv <- dtx[i,]	# actual deaths
      ev <- etx[i,]	# exposure
      kappa2v[i] <- llmaxM2D(beta1v, beta2v, beta3v,
                          kappa2v[i], gamma3v[(n+i-1):i], dv, ev, wv = wa[i,])}

    mhat <- mtx*0
    for(i in 1:m)
      mhat[i,] <- exp(beta1v + beta2v*kappa2v[i] + beta3v*gamma3v[(n+i-1):i])

    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat)
    l1      <- sum((dtx*log(etx*mhat) - etx*mhat - lgamma(dtx + 1))*wa)

    # Now apply the constraints
    if(constraints == "ALL"){
      fac21   <- mean(kappa2v)
      fac22   <- sqrt(sum(beta2v^2))
      kappa2v <- fac22*(kappa2v - fac21)     # sum kappa2=0
      beta2v  <- beta2v/fac22                # sum beta2=1
      beta1v  <- beta1v + beta2v*fac22*fac21 # => adjust beta1
    }
    else{
      fac22   <- sqrt(sum(beta2v^2))
      kappa2v <- fac22*kappa2v      # => adjust kappa2
      beta2v  <- beta2v/fac22  		  # sum beta2=1
    }

    mhat <- mtx*0
    for(i in 1:m)
      mhat[i,] <- exp(beta1v + beta2v*kappa2v[i] + beta3v*gamma3v[(n+i-1):i])

    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat)
    l1      <- sum((dtx*log(etx*mhat) - etx*mhat - lgamma(dtx + 1))*wa)

    # Stage 1A optimise over the beta1(x)
    if(all(int == 0)){
      for(j in 1:n){
        # cycle through the range of ages
        wv        <- 1	    # can be set to a vector of weights to e.g. exclude duff years
        wv        <- wa[,j]
        s1        <- sum(wv*dtx[,j])
        s2        <- sum(wv*etx[,j]*exp(beta2v[j]*kappa2v + beta3v[j]*gamma3v[(n+1-j):(n+m-j)]))
        beta1v[j] <- log(s1) - log(s2)
      }
    }

    mhat <- mtx*0
    for(i in 1:m)
      mhat[i,] <- exp(beta1v + beta2v*kappa2v[i] + beta3v*gamma3v[(n+i-1):i])

    epsilon <- (dtx - etx*mhat)/sqrt(etx*mhat)
    l1      <- sum((dtx*log(etx*mhat) - etx*mhat - lgamma(dtx + 1))*wa)
  }		 # end while loop

  # calculate number of parameters and deduct the number of constraints
  # also count number of parameters in 'beta1v'!
  if(constraints == "ALL")
    npar <- length(beta1v) + length(beta2v) + length(kappa2v) - 2 else
      npar <- length(beta1v) + length(beta2v) + length(kappa2v) - 1

  # Calculate the BIC
  BIC <- -2*l1 + log(sum(wa))*npar

  list(beta1 = beta1v, beta2 = beta2v, beta3 = beta3v,
       kappa2 = kappa2v, gamma3 = gamma3v, x = xv, y = yv, cy = cy,
       wa = wa, epsilon = epsilon, mhat = mhat, ll = l1, BIC = BIC, npar = npar,
       mtxLastYear = mtx[m,])
}


#' @keywords internal
llmaxM2B <- function(b1, b2, b3, k2, g3, dv, ev, wv = 1){
  #   b1,b3,k2,g3 are given
  #   solve for b2
  if(min(dv=0)) {
    print(b1)
    print(b2)
    print(b3)
    print(k2)
    print(g3)
    print(dv)
    print(ev)}

  b21    <- b2
  b20    <- b21-1
  thetat <- k2*ev*exp(b1 + b3*g3)
  s1     <- sum(dv*k2*wv)

  while(abs(b21 - b20) > 1e-10)
  {
    b20 <- b21;
    f0  <- sum((exp(b20*k2)*thetat)*wv) - s1;
    df0 <- sum((exp(b20*k2)*k2*thetat)*wv)
    b21 <- b20 - f0/df0
  }
  b21
}

#' @keywords internal
llmaxM2D <- function(b1, b2, b3, k2, g3, dv, ev, wv = 1){
  #   b1,b2,b3,g3 are given
  #   solve for k2
  k21    <- k2
  k20    <- k21-1
  thetat <- b2*ev*exp(b1 + b3*g3)
  s1     <- sum(dv*b2*wv)
  while(abs(k21 - k20) > 1e-10)
  {
    k20 <- k21
    f0  <- sum((exp(k20*b2)*thetat)*wv) - s1
    df0 <- sum((exp(k20*b2)*b2*thetat)*wv)
    k21 <- k20 - f0/df0
  }
  k21
}

