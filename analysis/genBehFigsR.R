

filesToLoad = 5:37

## What analyses are important?
# - Do errors during a dual trial cause the trial to fail more often
# - Do errors during a previous trial cause the current trial to fail
#
## Performance and Reaction time
#
# - Tasks, single, dual
# - Location
# - Gender
# - Order
# 
# - 

for fi = filesToLoad {
  # Dealing with run # fi
  main = read.csv(sprintf('../analysis/s300/csv/main%02.f.csv',fi))
  per = read.csv(sprintf('../analysis/s300/csv/per%02.f.csv',fi)
                 
  # Let's do some analysis!
  
  # First, let's clean this data.frame
  main = cleanup(main)
  per = cleanup(per)
}

cleanup <- function(data) {
  
}