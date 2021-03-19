#Create directories to store simulation results



if (!dir.exists("own_results")){
  succes<-dir.create("own_results")
  if (!succes){
    print("error in create directory")
  }
}


