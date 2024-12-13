
# Get number of samples by pop --------------------------------------------
ALL_SAMPLES <- "/project/lbernatchez/users/lalec31/projets_labo/Bastien/GTseq_saal_202410/dev_lcSNPs_for_GTseq/02_infos/samples_Bastien_Xavier_ID_POP_TYPE.txt"

all_samples <- read.table(ALL_SAMPLES, col.names = c('ID', 'POP', 'TYPE')) [, 1:2]

table(all_samples$POP)



# Split into a training and a test sets -----------------------------------
MIN_POP_SIZE <- 20  # min number of samples in a given pop for it to be sampled for a test subset
TO_SAMPLE <- 5      # number of samples to sample from each pop with sufficient MIN_POP_SIZE


test_set <- character()
training_set <- character()

for (pop in unique(all_samples$POP)){
  pop_set <- subset(all_samples, POP == pop)
  # Check number of samples available for pop
  if (nrow(pop_set) >= MIN_POP_SIZE){
    # Substract samples for testing
    test_IDs <- sample(x = pop_set$ID, size = TO_SAMPLE, replace = FALSE)
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
