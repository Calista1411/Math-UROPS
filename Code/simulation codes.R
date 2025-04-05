suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

f <- "all_positions.txt"
snp <- data.table(t(read.table(f, header=FALSE, fill=TRUE, row.names=NULL, sep = " "))) %>%
  rename(position = V1) %>%
  mutate(mut = 1)
filler <- data.frame(position = 1:29903)
snp <- right_join(snp, filler)
snp$mut[is.na(snp$mut)] <- 0

sum(snp$mut==1)

vcf <- fread("problematic_sites_sarsCov2.vcf", check.names=TRUE, skip="#CHROM")
mask <- c(vcf$POS)
# Additional sites to mask, from van Dorp et al. (2020, Supplementary Table S5):
# https://doi.org/10.1016/j.meegid.2020.104351
homoplasic <- c(11083,13402,21575,16887,27384,3037,25563,8782,10323,11074,14408,6255,12880,21137,1059,1820,1912,6040,15324,29353,14805,15720,18060,28077,28826,28887,29253,29553,29700,6310,6312,11704,14786,17747,18756,20148,22661,23010,23403,29095,29422,3177,4084,6990,8078,11916,14724,14925,17247,18788,18877,20755,21648,24034,25947,26152,26461,27005,27046,27964,28881,29742,1457,4255,5784,7011,8293,8917,9223,10319,10507,11320,12781,13947,15760,16260,19684,22988,23422,24390,25916,26144,26530,26730,27525,28144,28311,28344,28851,28854,29751,379,541,833,884,1076,1570,1594,2113,3096,3253,3787,4113,4320,6573,7438,7765,10789,10851,11417,14747,15960,16762,17410,17639,17799,17858,18656,20031,20268,20275,21204,21707,23533,23587,24368,24389,24694,25494,25688,26211,26729,26735,28657,28688,28739,28857,28878,29540,29585,29734,313,490,1515,2455,2558,4809,6723,7479,8767,9477,9479,10097,10265,10450,11195,11801,13730,13929,14741,14912,15277,15927,16289,16381,17104,17373,17690,17944,18652,18713,18928,18998,19170,20931,23086,23707,23731,23929,24054,24862,25433,25572,25979,26124,26625,26936,27299,27635,27679,28580,28821,28836,28882,28883,29144,29635,29686)
mask <- c(mask, homoplasic)
mask <- unique(mask)

# Add bases from reference sequence
g <- data.frame(base = t(read.table("refseq.txt", header=FALSE, fill=TRUE, row.names=NULL, sep = " ", colClasses="character"))) %>%
  mutate(position = 1:29904)
snp <- right_join(snp, g)
snp <- snp[-29904,]

# Mask the sites
snp$mut[mask] <- NA

suppressWarnings(suppressMessages(library(wavethresh)))
suppressWarnings(suppressMessages(library(EbayesThresh)))

# Carry out the wavelet regression, splitting by base type
mut_probs <- data.frame(position = 1:29903, rate = rep(NA, 29903), base = as.character(snp$base), stringsAsFactors=FALSE)
snp$mut[is.na(snp$mut)] <- 0
base_to_num <- c(A = 1, C = 2, G = 3, T = 4)
# Iterate through each base type
for(b in c("A", "C", "G", "T")) {
  # Filter the data for the current base
  X <- snp$base[snp$base == b]
  Xpos <- snp$position[snp$base == b]
  Xpos <- Xpos[!is.na(X)]  # Remove NA positions
  X <- X[!is.na(X)]  # Remove NA values for bases
  # Convert bases to numeric values for processing
  X_numeric <- base_to_num[X]
  # Padding to make total length a power of 2
  pow <- ceiling(log(length(X_numeric), base = 2))  
  desired_length <- 2^pow  
  # Reflecting the data at the endpoints
  ext_left <- floor((desired_length - length(X_numeric)) / 2)
  ext_right <- desired_length - length(X_numeric) - ext_left
  # Adjust the padding, ensuring it doesn't go out of bounds
  Xf <- c(rev(X_numeric[1:ext_left]), X_numeric, rev(X_numeric[(length(X_numeric) - ext_right + 1):length(X_numeric)]))
  # Check that the length of Xf is a power of 2
  print(length(Xf))  # This should now be a power of 2
  
  # Wavelets with empirical Bayes thresholding
  Y <- wd(Xf, filter.number=6, family="DaubLeAsymm")
  print(str(Y))
  print(sum(is.na(Y$D)))
  #Y <- ebayesthresh.wavelet(Y, threshrule="soft")
  Y$D <- sapply(Y$D, function(x) ebayesthresh(x, threshrule="soft"))
  Y <- wr(Y)
  Y <- Y[(ext_left+1):(ext_left + length(X_numeric))]
  # Check if Y has valid values
  if (length(Y) == 0 || all(is.na(Y))) {
    print(paste("No valid Y values for base", b))
  } else {
    print(paste("Valid Y values for base", b))
    
    # Update mut_probs$rate for positions corresponding to the base
    mut_probs$rate[mut_probs$position %in% Xpos] <- Y
  }
}

mut_probs <- mut_probs %>%
  mutate(rate = ifelse(rate < 0, 0, rate),
         rate_norm = rate/sum(rate, na.rm=TRUE))

# Define color palette
color_palette <- c("A" = "#F7C6C7", "C" = "#B6E2D3", "G" = "#A8C8F0", "T" = "#D6C6E1")
color_palette2 <- c(
  A = "#E88C8D",  # Deeper blush pink
  C = "#78C6B0",  # Rich mint green
  G = "#6FA4DC",  # Softer sky blue with more saturation
  T = "#A38BC6"   # Medium lavender
)

library(ggplot2)
library(latex2exp)
library(scales)

# Set plot dimensions 
options(repr.plot.width = 15, repr.plot.height = 6)

# Fonts
library(extrafont)
library(showtext)
font_add(family = "CM Roman", regular = "/Users/devinacalista/Library/Fonts/cmunrm copy.ttf")  # Use the correct path
showtext_auto()

# Graph of mutation rates density across the genome 
g_smooth <- ggplot(data = mut_probs, aes(x = position, y = rate_norm)) +
  geom_smooth(aes(colour = base), method = "loess", span = 0.05, size = 1, se = FALSE) +
  scale_color_manual(values = color_palette2) +  
  theme_minimal() +
  xlab("Position") +
  ylab(expression(tilde(P))) +
  scale_x_continuous(
    breaks = seq(0, 30000, 5000),
    labels = seq(0, 30000, 5000),
    minor_breaks = seq(0, 30000, 1000)
  ) +
  scale_y_continuous(
    breaks = seq(1.317471e-05, 5.269884e-05, length.out = 5),
    labels = scales::label_number(scale = 1, accuracy = 1e-5),
    minor_breaks = NULL
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "CM Roman")
  ) +
  ggtitle("Smoothed Mutation Rates by Position")

g_smooth

# Algorithm 1
sim <- function(M, m, P) {
  unique_sites <- integer(m)   # Store unique sites
  seen_sites <- integer(M)     # Track seen sites 
  unique_count <- 0            # Counter for unique sites
  step <- 0                    # Total number of runs or steps taken
  
  while (unique_count < m && step < 100000) {  # Prevent infinite loop
    step <- step + 1
    visited <- sample(1:M, size = 1, replace = TRUE, prob = P)
    
    if (seen_sites[visited] == 0) {  # New site found
      unique_count <- unique_count + 1
      unique_sites[unique_count] <- visited
      seen_sites[visited] <- 1
    }
  }
  
  return(step)  # Return the total number of steps taken
}

calc <- function(M, m, n, P) {
  results <- rep(NA, n)
  for(i in 1:n) {
    results[i] <- sim(M, m, P)
  }
  return(results - m)
}

Y <- mut_probs$rate_norm
m <- 861*1.1 # Number of variable sites/SNPs (obtained from snp-sites) * factor of 1.1
# 349 for Asia sample, 521 for Europe sample, 388 for Africa sample, 861 for North America sample
res <- calc(length(Y), m, 1000, Y/sum(Y))

# Plot histogram of the distribution of the number of recurrent mutations
options(repr.plot.width=10, repr.plot.height=5)
h1 <- ggplot(data.frame(x = res), aes(x = x)) +
  geom_histogram(aes(y = ..count.. / sum(..count..)), 
                 binwidth = 1, colour = "purple", fill = "pink", alpha = 0.5) +
  geom_text(stat = "bin",
            aes(y = ..count.. / sum(..count..), 
                label = ifelse(..count.. > 0, 
                               scales::percent(..count.. / sum(..count..), accuracy = 0.1), "")),
            binwidth = 1,
            vjust = -1,
            size = 2) +
  theme_minimal() +
  xlab("Number of recurrent mutations") +
  ylab("Proportion") +
  scale_x_continuous(breaks = seq(0, 40, 5), labels = seq(0, 40, 5)) +
  ggtitle("Simulated Distribution of Number of Recurrent Mutations") +
  theme(text = element_text(family = "CM Roman"))

h1