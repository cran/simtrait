## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- cache = FALSE, include = FALSE------------------------------------------
## copied from examples from the "simmer" R package
## after: https://www.enchufa2.es/archives/suggests-and-vignettes.html
## by Iñaki Úcar
required <- c("popkin", "bnpsd") # only suggested since simtrait doesn't need them to run...

if (!all(sapply(required, requireNamespace, quietly = TRUE)))
  knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
library(popkin)   # to create plots of our covariance matrices
library(bnpsd)    # to simulate an admixed population
library(simtrait) # this package

## -----------------------------------------------------------------------------
# dimensions of data/model
# number of loci
m_loci <- 10000
# number of individuals, smaller than usual for easier visualizations
n_ind <- 30
# number of intermediate subpops (source populations for admixed individuals)
k_subpops <- 3

# define population structure
# FST values for 3 subpopulations (proportional/unnormalized)
inbr_subpops <- 1 : k_subpops
bias_coeff <- 0.5 # bias coeff of standard Fst estimator
Fst <- 0.3 # desired final Fst
obj <- admix_prop_1d_linear(
    n_ind = n_ind,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    coanc_subpops = inbr_subpops,
    fst = Fst
)
admix_proportions <- obj$admix_proportions
# rescaled Fst vector for intermediate subpops
inbr_subpops <- obj$coanc_subpops

# get pop structure parameters of the admixed individuals
coancestry <- coanc_admix(admix_proportions, inbr_subpops)
kinship <- coanc_to_kinship(coancestry)

# draw allele freqs and genotypes
out <- draw_all_admix(admix_proportions, inbr_subpops, m_loci)
X <- out$X # genotypes
p_anc <- out$p_anc # ancestral AFs

## -----------------------------------------------------------------------------
# parameters of simulation
m_causal <- 100
herit <- 0.8
# default 0, let's try a non-trivial case
mu <- 1
# default 1, also let's see that this more complicated case works well
sigma_sq <- 1.5

# create simulated trait
# case of exact p_anc
obj <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
    p_anc = p_anc,
    mu = mu,
    sigma_sq = sigma_sq
)
# trait vector
length(obj$trait)
n_ind
obj$trait
# randomly-picked causal locus indexes
length( obj$causal_indexes )
m_causal
head( obj$causal_indexes ) # show partially...
# regression coefficients vector
length( obj$causal_coeffs )
m_causal
head( obj$causal_coeffs ) # show partially...

## -----------------------------------------------------------------------------
# the theoretical covariance matrix of the trait is calculated by cov_trait
V <- cov_trait(kinship = kinship, herit = herit, sigma_sq = sigma_sq)

# simulate these many traits
n_traits <- 1000
# store in this matrix, initialize with zeroes
Y_rc_freq <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
# start loop
for (i in 1 : n_traits) {
    obj <- sim_trait(
        X = X,
        m_causal = m_causal,
        herit = herit,
        p_anc = p_anc,
        mu = mu,
        sigma_sq = sigma_sq
    )
    Y_rc_freq[i,] <- obj$trait # store in i^th row
}
# estimate sample covariance
V_rc_freq <- cov(Y_rc_freq)

## ---- fig.width = 3, fig.align = 'center'-------------------------------------
par_orig <- par(mgp = c(2, 0.5, 0))
# reduce margins from default
par(mar = c(3.5, 3, 0, 0) + 0.2)
# visualize distribution
boxplot(
    list(
        'RC freq' = rowMeans(Y_rc_freq)
    ),
    xlab = "Trait Type",
    ylab = 'Sample Mean'
)
# red line marks expected mean
abline(h = mu, col = 'red')
par( par_orig ) # reset `par`

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
plot_popkin(
    list(V, V_rc_freq),
    titles = c('Theoretical', 'RC freq'),
    leg_title = 'Covariance',
    panel_letters_adj = 0,
    # set margin for title (top is non-zero)
    mar = c(0, 2)
)

## -----------------------------------------------------------------------------
# store this in new matrix
Y_rc_kin <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
# start loop
for (i in 1 : n_traits) {
    obj <- sim_trait(
        X = X,
        m_causal = m_causal,
        herit = herit,
        # whole kinship matrix can be passed instead of just mean
        kinship = kinship,
        mu = mu,
        sigma_sq = sigma_sq
    )
    Y_rc_kin[i,] <- obj$trait # store in i^th row
}
# estimate sample covariance
V_rc_kin <- cov(Y_rc_kin)

## ---- fig.width = 3, fig.align = 'center'-------------------------------------
par_orig <- par(mgp = c(2, 0.5, 0))
# reduce margins from default
par(mar = c(3.5, 3, 0, 0) + 0.2)
# visualize distribution
boxplot(
    list(
        "RC freq" = rowMeans(Y_rc_freq),
        "RC kinship" = rowMeans(Y_rc_kin)
    ),
    xlab = "Trait Type",
    ylab = 'Sample Mean'
)
# red line marks expected mean
abline(h = mu, col = 'red')
par( par_orig ) # reset `par`

## ---- fig.width = 7, fig.height = 2.35, fig.align = 'center'------------------
plot_popkin(
    list(V, V_rc_freq, V_rc_kin),
    titles = c('Theoretical', 'RC freq', 'RC kinship'),
    leg_title = 'Covariance',
    panel_letters_adj = 0,
    mar = c(0, 2)
)

## -----------------------------------------------------------------------------
# This function simulates trait replicates in one call,
# generating a matrix comparable to the previous ones.
Y_mvn <- sim_trait_mvn(
    rep = n_traits,
    kinship = kinship,
    herit = herit,
    mu = mu,
    sigma_sq = sigma_sq
)
# estimate sample covariance
V_mvn <- cov(Y_mvn)

## ---- fig.width = 4, fig.align = 'center'-------------------------------------
par_orig <- par(mgp = c(2, 0.5, 0))
# reduce margins from default
par(mar = c(3.5, 3, 0, 0) + 0.2)
# visualize distribution
boxplot(
    list(
        "RC freq" = rowMeans(Y_rc_freq),
        "RC kinship" = rowMeans(Y_rc_kin),
        "MVN" = rowMeans(Y_mvn)
    ),
    xlab = "Trait Type",
    ylab = 'Sample Mean'
)
# red line marks expected mean
abline(h = mu, col = 'red')
par( par_orig ) # reset `par`

## ---- fig.width = 7, fig.height = 2, fig.align = 'center'---------------------
plot_popkin(
    list(V, V_rc_freq, V_rc_kin, V_mvn),
    titles = c('Theoretical', 'RC freq', 'RC kinship', 'MVN'),
    leg_title = 'Covariance',
    panel_letters_adj = 0,
    mar = c(0, 2),
    leg_width = 0.4
)

## -----------------------------------------------------------------------------
# store this in new matrix
Y_fes_freq <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
Y_fes_kin <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
# start loop
for (i in 1 : n_traits) {
    obj <- sim_trait(
        X = X,
        m_causal = m_causal,
        herit = herit,
        p_anc = p_anc,
        mu = mu,
        sigma_sq = sigma_sq,
        fes = TRUE # only diff from orig run
    )
    Y_fes_freq[i,] <- obj$trait # store in i^th row

    obj <- sim_trait(
        X = X,
        m_causal = m_causal,
        herit = herit,
        kinship = kinship,
        mu = mu,
        sigma_sq = sigma_sq,
        fes = TRUE # only diff from orig run
    )
    Y_fes_kin[i,] <- obj$trait # store in i^th row
}
# estimate sample covariance
V_fes_freq <- cov(Y_fes_freq)
V_fes_kin <- cov(Y_fes_kin)

## ---- fig.width = 6, fig.align = 'center'-------------------------------------
par_orig <- par(mgp = c(2, 0.5, 0))
# reduce margins from default
par(mar = c(3.5, 3, 0, 0) + 0.2)
# visualize distribution
boxplot(
    list(
        "RC freq" = rowMeans(Y_rc_freq),
        "RC kinship" = rowMeans(Y_rc_kin),
        "MVN" = rowMeans(Y_mvn),
        "FES freq" = rowMeans(Y_fes_freq),
        "FES kinship" = rowMeans(Y_fes_kin)
    ),
    xlab = "Trait Type",
    ylab = 'Sample Mean'
)
# red line marks expected mean
abline(h = mu, col = 'red')
par( par_orig ) # reset `par`

## ---- fig.width = 7, fig.height = 4, fig.align = 'center'---------------------
plot_popkin(
    list( V, V_rc_freq, V_rc_kin, V_mvn, V_fes_freq, V_fes_kin ),
    titles = c('Theoretical', 'RC freq', 'RC kinship', 'MVN', 'FES freq', 'FES kinship'),
    leg_title = 'Covariance',
    panel_letters_adj = 0,
    mar = c(0, 2),
    leg_width = 0.4,
    layout_rows = 2
)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'---------------------
par_orig <- par(mgp = c(2, 0.5, 0))
# create multipanel figure
par( mfrow = c(2, 3) )
# reduce margins from default
par(mar = c(3.5, 3, 0, 0) + 0.2)
plot( V, V_rc_freq, xlab = 'Theoretical Cov', ylab = 'RC freq Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_rc_kin, xlab = 'Theoretical Cov', ylab = 'RC kinship Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_mvn, xlab = 'Theoretical Cov', ylab = 'MVN Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_fes_freq, xlab = 'Theoretical Cov', ylab = 'FES freq Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_fes_kin, xlab = 'Theoretical Cov', ylab = 'FES kinship Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
barplot(
    c(
        rmsd( V, V_rc_freq ),
        rmsd( V, V_rc_kin ),
        rmsd( V, V_mvn ),
        rmsd( V, V_fes_freq ),
        rmsd( V, V_fes_kin )
    ),
    names.arg = c('RC freq', 'RC kin', 'MVN', 'FES freq', 'FES kin'),
    ylab = 'RMSD from Theoretical',
    las = 3
)
par( par_orig ) # reset `par`

## -----------------------------------------------------------------------------
# first level, first half is all one group called "a"
# second half is another group called "b"
n_half <- round( n_ind / 2 )
labs1 <- c(
    rep.int( 'a', n_half ),
    rep.int( 'b', n_ind - n_half )
)
# second level will be in thirds instead
# each level is independent, so group names can repeat across levels
n_third <- round( n_ind / 3 )
labs2 <- c(
    rep.int( 'a', n_third ),
    rep.int( 'b', n_third ),
    rep.int( 'c', n_ind - 2 * n_third )
)
# combine!
labs <- cbind( labs1, labs2 )
# set heritability and variance effects for these levels
# reduce heritability to allow for large group variances
herit <- 0.3
labs_sigma_sq <- c( 0.3, 0.2 )

## -----------------------------------------------------------------------------
V <- cov_trait( kinship = kinship, herit = herit, sigma_sq = sigma_sq,
               labs = labs, labs_sigma_sq = labs_sigma_sq )
# store in this matrix, initialize with zeroes
Y_rc_freq  <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
Y_rc_kin   <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
Y_fes_freq <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
Y_fes_kin  <- matrix(data = 0, nrow = n_traits, ncol = n_ind)
# start loop
for (i in 1 : n_traits) {
    Y_rc_freq[i,]  <- sim_trait( X=X, m_causal=m_causal, herit=herit, mu=mu, sigma_sq=sigma_sq, 
        labs=labs, labs_sigma_sq=labs_sigma_sq, p_anc=p_anc )$trait
    Y_rc_kin[i,]   <- sim_trait( X=X, m_causal=m_causal, herit=herit, mu=mu, sigma_sq=sigma_sq, 
        labs=labs, labs_sigma_sq=labs_sigma_sq, kinship=kinship )$trait
    Y_fes_freq[i,] <- sim_trait( X=X, m_causal=m_causal, herit=herit, mu=mu, sigma_sq=sigma_sq, 
        labs=labs, labs_sigma_sq=labs_sigma_sq, p_anc=p_anc, fes=TRUE )$trait
    Y_fes_kin[i,]  <- sim_trait( X=X, m_causal=m_causal, herit=herit, mu=mu, sigma_sq=sigma_sq, 
        labs=labs, labs_sigma_sq=labs_sigma_sq, kinship=kinship, fes=TRUE )$trait
}
Y_mvn <- sim_trait_mvn( rep = n_traits, kinship = kinship, herit = herit,
    mu = mu, sigma_sq = sigma_sq, labs = labs, labs_sigma_sq = labs_sigma_sq )
# estimate sample covariance
V_rc_freq  <- cov(Y_rc_freq)
V_rc_kin   <- cov(Y_rc_kin)
V_fes_freq <- cov(Y_fes_freq)
V_fes_kin  <- cov(Y_fes_kin)
V_mvn      <- cov(Y_mvn)

## ---- fig.width = 6, fig.align = 'center'-------------------------------------
par_orig <- par(mgp = c(2, 0.5, 0))
# reduce margins from default
par(mar = c(3.5, 3, 0, 0) + 0.2)
# visualize distribution
boxplot(
    list(
        "RC freq" = rowMeans(Y_rc_freq),
        "RC kinship" = rowMeans(Y_rc_kin),
        "MVN" = rowMeans(Y_mvn),
        "FES freq" = rowMeans(Y_fes_freq),
        "FES kinship" = rowMeans(Y_fes_kin)
    ),
    xlab = "Trait Type",
    ylab = 'Sample Mean'
)
# red line marks expected mean
abline(h = mu, col = 'red')
par( par_orig ) # reset `par`

## ---- fig.width = 7, fig.height = 4, fig.align = 'center'---------------------
plot_popkin(
    list( V, V_rc_freq, V_rc_kin, V_mvn, V_fes_freq, V_fes_kin ),
    titles = c('Theoretical', 'RC freq', 'RC kinship', 'MVN', 'FES freq', 'FES kinship'),
    leg_title = 'Covariance',
    panel_letters_adj = 0,
    mar = c(0, 2),
    leg_width = 0.4,
    layout_rows = 2
)

## ---- fig.width = 7, fig.height = 5, fig.align = 'center'---------------------
par_orig <- par(mgp = c(2, 0.5, 0))
# create multipanel figure
par( mfrow = c(2, 3) )
# reduce margins from default
par(mar = c(3.5, 3, 0, 0) + 0.2)
plot( V, V_rc_freq, xlab = 'Theoretical Cov', ylab = 'RC freq Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_rc_kin, xlab = 'Theoretical Cov', ylab = 'RC kinship Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_mvn, xlab = 'Theoretical Cov', ylab = 'MVN Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_fes_freq, xlab = 'Theoretical Cov', ylab = 'FES freq Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
plot( V, V_fes_kin, xlab = 'Theoretical Cov', ylab = 'FES kinship Cov' )
abline( 0, 1, lty = 2, col = 'gray' )
barplot(
    c(
        rmsd( V, V_rc_freq ),
        rmsd( V, V_rc_kin ),
        rmsd( V, V_mvn ),
        rmsd( V, V_fes_freq ),
        rmsd( V, V_fes_kin )
    ),
    names.arg = c('RC freq', 'RC kin', 'MVN', 'FES freq', 'FES kin'),
    ylab = 'RMSD from Theoretical',
    las = 3
)
par( par_orig ) # reset `par`

