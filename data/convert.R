source("schools.data.R")
library(jsonlite)

data <- list(
  Gender = Gender,
  LRT = LRT,
  M = M,
  N = N,
  R = R,
  school = school,
  School_denom = School_denom,
  School_gender = School_gender,
  Y = Y,
  VR = VR
)

json <- toJSON(data, pretty = TRUE)
output_filename <- "data.json"
cat(json, file = output_filename, append = TRUE)