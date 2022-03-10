library(jsonlite)

list2json <- function(data, path) {
  json <- toJSON(data, pretty = TRUE, auto_unbox = TRUE) # unbox vector of length 1
  cat(json, file = path, append = TRUE)
}

# Data
source("schools.data.R")
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
list2json(data, "data.json")

# Init

## Init 1
source("init1.R")
list2json(init, "init1.json")

## Init 2
source("init2.R")
list2json(init, "init2.json")