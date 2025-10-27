.onLoad <- function(...){
  packageStartupMessage("\n")
  packageStartupMessage("Welcome to SpLin.")
  packageStartupMessage("\n")
  packageStartupMessage("Version: ",utils::packageDescription('SpLin')$Version)
  packageStartupMessage("\n")
  packageStartupMessage("If this is your first time running SpLin you should see ?SpLin")
  packageStartupMessage("\n")
}