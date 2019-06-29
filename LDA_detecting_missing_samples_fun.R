LDA_detecting_missing_samples_fun <- function(Sample_name_endings,colnames_data2m,nb_repeats,byRow){

NEXT_CHAR <- function(previous_last_char_index) {
  if (previous_last_char_index == length(Sample_name_endings)) 
  {
    next_last_char_index = 1
  } else
    next_last_char_index = previous_last_char_index + 1
  next_last_char_index
}
# Missing_positions_indices is a vector of the same length as the number of samples (missing or not) in the grid,
# taking value 1 for a missing sample and 0 otherwise
Missing_positions_indices = vector(length=length(Sample_name_endings)*nb_repeats,mode="numeric")
# Dealing with the case of missing samples at the beginning of the grid:
char_string = colnames_data2m[1]
last_char = substr(char_string,nchar(char_string),nchar(char_string))
last_char_index = which(Sample_name_endings==last_char)
if (last_char_index!=1)
{
  spatial_index = last_char_index
  Missing_positions_indices[1:(last_char_index-1)] = 1
} else
  spatial_index = 1
previous_last_char_index = spatial_index  
# Loop over non-empty samples, i.e. over the columns of data2m:
for (j in 2:length(colnames_data2m))
{
  spatial_index = spatial_index+1
  char_string = colnames_data2m[j]
  last_char = substr(char_string,nchar(char_string),nchar(char_string))
  last_char_index = which(Sample_name_endings==last_char)
  expected_char = NEXT_CHAR(previous_last_char_index)
  # If the observed name ending is not the expected one,
  # Missing_positions_indices is filled with as many 1 as there are missing samples:
  while (last_char_index != expected_char)
  {
    Missing_positions_indices[spatial_index] = 1
    expected_char = NEXT_CHAR(expected_char)
    spatial_index = spatial_index+1
  }
  previous_last_char_index = last_char_index
}
# Dealing with the case of missing samples at the end of the grid:
if (previous_last_char_index != length(Sample_name_endings))
{
  last_char_index = length(Sample_name_endings)
  expected_char = previous_last_char_index
  while (last_char_index != expected_char)
  {
    spatial_index = spatial_index+1
    Missing_positions_indices[spatial_index] = 1
    expected_char = NEXT_CHAR(expected_char)
  }
}

return(Missing_positions_indices)
}