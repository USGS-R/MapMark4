current_wd <- setwd("data-raw")

tmp <- read.csv( "ExampleGatm.csv", header = TRUE)

ExampleTonnageData <-
  data.frame(Name = tmp$Deposit.Name,
             Ore = tmp$Tonnage..mt.,
             Cu = tmp$Tonnage..mt. * tmp$Cu.Grade.... * 0.01,
             Au = tmp$Tonnage..mt. * tmp$Au.Grade.... * 0.01)

devtools::use_data(ExampleTonnageData, overwrite = TRUE)

setwd(current_wd)
