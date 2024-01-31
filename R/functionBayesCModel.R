## Function for BayesC model 

buildBayesCModel <- function(trainGenoTable,trainPhenoTable,testGenoTable,testPhenoTable,heritability){ 

	
	 trainPheno <- trainPhenoTable[,1]
         trainGeno <- trainGenoTable
	 h2 <- heritability
	 nIndividuals <- 2000
	 
	 
	predBayesC <- BGLR(trainPheno,response_type="gaussian",ETA=list(list(X=trainGeno,model="BayesC")),nIter=41000,burnIn=1000,df0=4,R2=h2)
	
	
###############################################################################3

#################################################################################################################
  
  testPheno <- testPhenoTable[,1] 
  testGeno <- testGenoTable 
  
  train_nIndividuals <- (0.8*nIndividuals)
  test_nIndividuals <-  (0.2*nIndividuals)
  
  Pred.pheno.train <- (trainGeno %*% predBayesC$ETA[[1]]$b)
  Pred.pheno.valid <- (testGeno %*% predBayesC$ETA[[1]]$b) 
  
  #Pred.pheno.valid <- predBL_BRR$mu+(as.vector(pred.valid))
  
  Pred.Accuracy <- cor(Pred.pheno.valid,testPheno)
  
  train_MSE <- (sum((Pred.pheno.train-trainPheno))^2)/(train_nIndividuals)
  test_MSE <- (sum((Pred.pheno.valid-testPheno))^2)/(test_nIndividuals) 
  
  sol <- list(Pred.Accuracy,predBayesC,predBayesC$ETA[[1]]$b,train_MSE,test_MSE)
  
  return(sol)


}
