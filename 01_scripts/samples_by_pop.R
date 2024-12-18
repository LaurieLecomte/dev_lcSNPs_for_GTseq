
# Get number of samples by pop --------------------------------------------

argv <- commandArgs(T)
ALL_SAMPLES <- argv[1]
#ALL_SAMPLES <- "02_infos/ID_POP.txt"

all_samples <- read.table(ALL_SAMPLES, col.names = c('ID', 'POP'))

table(all_samples$POP)

set.seed(42)

# Split into a training and a test sets -----------------------------------
MIN_POP_SIZE <- 20  # min number of samples in a given pop for it to be sampled for a test subset
TO_SAMPLE <- 5      # number of samples to sample from each pop with sufficient MIN_POP_SIZE


test_set <- character()
training_set <- character()

for (pop in unique(all_samples$POP)){
  pop_set <- subset(all_samples, POP == pop)
  pop_size <- nrow(pop_set)
  to_sample <- min(MIN_SAMPLES, (pop_size - MIN_POP_SIZE))
  
  # Check number of samples available for pop
  if (nrow(pop_set) >= MIN_POP_SIZE){
    # Substract samples for testing
    test_IDs <- sample(x = pop_set$ID, size = to_sample, replace = FALSE)
    training_IDs <- pop_set$ID[!pop_set$ID %in% test_IDs]
    
    test_set <- append(test_set, test_IDs)
    training_set <- append(training_set, training_IDs)
  } else {
    # Use all avail samples as training set if not enough samples avail for pop
    training_IDs <- pop_set$ID
    training_set <- append(training_set, training_IDs)
  }
}

# Export to lists
writeLines(text = test_set, con = '02_infos/test_samples.txt')
writeLines(text = training_set, con = '02_infos/training_samples.txt')
