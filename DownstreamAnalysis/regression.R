
library(ggplot2)
library(reshape2)
library(glmnet)
#################
## normaliz the data
###############
standardize <- function(x) {
      if (!is.matrix(x)) dim(x) <- c(length(x),1)
        means <- colMeans(x)
        x <- t(x) - means
        vars <- sqrt(rowSums(x^2))
        x <- t(x/vars)
        #return(list(x=x, center=means, scale=vars))
        return(x)
      }
##################################
### Regression
##################################
LassoRegression = function(Xtrain, Y, alphaGiven,numberIter=1 , Output, verbose = TRUE,sampleName){
      # Y=(Y-mean(Y))/sd(Y)
      Xtrain=standardize(Xtrain)
      Xtrain[is.na(Xtrain)]=0
      sampleName=rownames(Xtrain)
      sampleSize=length(Y)

      predictions=rep(0,sampleSize)
      labs={}
      RSquare={}

      NonZeroGenes = {}
      NonZeroCoefficient = {}

      #Find best lambda
      cv.fit=cv.glmnet(Xtrain, Y, family="gaussian", type.measure = "deviance", nfolds = sampleSize,standardize = F)
      #plot(cv.fit)
      LambdaMin= cv.fit$lambda.min
      lambdas = cv.fit$lambda
      best.lam.ind = which(lambdas == LambdaMin)

      model=glmnet(Xtrain, Y, family="gaussian", lambda = lambdas , alpha=alphaGiven,standardize = F)

      for(i in 1:sampleSize)

            {




                    y <- predicted = c( predict(model, newx=as.matrix(t(Xtrain[i,])), s=LambdaMin ,type = "link"))

                    predictions[i]=y <- predicted
                    # Sum of Squares Total and Error

                    sse <- (y <- predicted - Y[i])^2
                    RSquare=c(RSquare,sse)

                  }


      NonZeroGenes = c(NonZeroGenes, colnames(Xtrain[,which(model$beta[, best.lam.ind] != 0)]))
      NonZeroCoefficient = c(NonZeroCoefficient, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])


      fileName=paste0("NonZeroGenes_Alpha_",alphaGiven,"_Sample_",sampleName[i],"_out")

      path=paste0(Output,alphaGiven,"/")
      dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

      #jpeg(file=paste(path,fileName,"_.jpg"),width = 780, height = 780,quality =75)
      #plotres(model)
      #dev.off()
      write.csv(NonZeroGenes, file=paste0(path,fileName,"_geneName.txt"))
      write.csv(NonZeroCoefficient, file=paste0(path,fileName,"_geneCoef.txt"))


      sst <- sum((Y - mean(Y))^2)
      # R squared
      sse=sum(RSquare)
      rsq <- 1 - sse / sst
      write.csv(rsq, file=paste0(path,"_RSQ.txt"))
      L=list(pred=predictions,rsq=rsq)
      class(L) = c("leave-oneOut")
      return(L)
    }
########################################
### Linear Regression Leave-one out function
########################################
LeaveOneOut = function(Xtrain, Y, alphaGiven,numberIter=1 , Output, verbose = TRUE,sampleName){
      # Y=(Y-mean(Y))/sd(Y)
      Xtrain=standardize(Xtrain)
      Xtrain[is.na(Xtrain)]=0
      sampleName=rownames(Xtrain)
      sampleSize=length(Y)

      predictions=rep(0,sampleSize)
      labs={}
      RSE={}
      RSTT={}
      for(i in 1:sampleSize)

            {


                    NonZeroGenes = {}
                    NonZeroCoefficient = {}

                    #Find best lambda
                    cv.fit=cv.glmnet(Xtrain[-i,], Y[-i], family="gaussian", type.measure = "deviance", nfolds = sampleSize-1,standardize=F)

                    LambdaMin= cv.fit$lambda.min
                    lambdas = cv.fit$lambda
                    best.lam.ind = which(lambdas == LambdaMin)

                    model=glmnet(Xtrain[-i,], Y[-i], family="gaussian", lambda = lambdas , alpha=alphaGiven,standardize=F)



                    y <- predicted = c( predict(model, newx=as.matrix(t(Xtrain[i,])), s=LambdaMin ,type = "link"))

                    predictions[i]=y <- predicted
                    # Sum of Squares Total and Error

                    sse <- (y <- predicted - Y[i])^2
                    RSE=c(RSE,sse)
                    RSTT=c(RSTT,sum((Y[-i]-mean(Y[-i]))^2))

                    NonZeroGenes = c(NonZeroGenes, colnames(Xtrain[,which(model$beta[, best.lam.ind] != 0)]))
                    NonZeroCoefficient = c(NonZeroCoefficient, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])

                    fileName=paste0("NonZeroGenes_Alpha_",alphaGiven,"_Sample_",sampleName[i],"_out")

                    path=paste0(Output,alphaGiven,"/")
                    dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

                    jpeg(file=paste(path,fileName,"_.jpg"),width = 780, height = 780,quality =75)
                    plot(cv.fit)
                    dev.off()

                    #jpeg(file=paste(path,fileName,"_.jpg"),width = 780, height = 780,quality =75)
                    #plotres(model)
                    #dev.off()
                    write.csv(NonZeroGenes, file=paste0(path,fileName,"_geneName.txt"))
                    write.csv(NonZeroCoefficient, file=paste0(path,fileName,"_geneCoef.txt"))


                  }


      # R squared
      sse=sum(RSE)
      rsq <- 1 - sse / mean(RSTT)
      write.csv(rsq, file=paste0(path,"_RSQ.txt"))
      L=list(pred=predictions,rsq=rsq)
      class(L) = c("leave-oneOut")
      return(L)
    }
########################################
### Linear Regression Leave-one out function
########################################
BootstrapwithCV = function(Xtrain, Y, alphaGiven,numberIter=1 , Output, verbose = TRUE,sampleName){
      sampleSize=length(Y)
        tmpX=Xtrain
        predictions=rep(0,sampleSize)
        labs={}
        RSquare={}
        sampleSize=length(Y)
        for(i in 1:1000)

              {

                      bootstrapIndex=sample(1:sampleSize,sampleSize, replace=T)
                      while(sd(bootstrapIndex)==0){
                                bootstrapIndex=sample(1:sampleSize,sampleSize, replace=T)
                                    }
                      Xtrain=standardize(tmpX[bootstrapIndex,])
                      Xtrain[is.na(Xtrain)]=0
                      sampleName=rownames(Xtrain)

                      cat(bootstrapIndex,"\n")
                      NonZeroGenes = {}
                      NonZeroCoefficient = {}

                      #Find best lambda
                      tryCatch({ cv.fit=cv.glmnet(Xtrain, Y[bootstrapIndex], family="gaussian", type.measure = "deviance", nfolds = sampleSize,standardize = F)
                              LambdaMin= cv.fit$lambda.min
                                  lambdas = cv.fit$lambda
                                  best.lam.ind = which(lambdas == LambdaMin)

                                  model=glmnet(Xtrain, Y[bootstrapIndex], family="gaussian", lambda = lambdas , alpha=alphaGiven,standardize = F)



                                  NonZeroGenes = c(NonZeroGenes, colnames(Xtrain[,which(model$beta[, best.lam.ind] != 0)]))
                                  NonZeroCoefficient = c(NonZeroCoefficient, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])

                                  fileName=paste0("NonZeroGenes_Alpha_",alphaGiven,"_Iteration_",i)

                                  path=paste0(Output,alphaGiven,"/")
                                  dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

                                  write.csv(NonZeroGenes, file=paste0(path,fileName,"_geneName.txt"))
                                  write.csv(NonZeroCoefficient, file=paste0(path,fileName,"_geneCoef.txt"))


                                  },error=function(e){})


                    }

      }
############################################################################################
# main
###############################################################################3

path="/H4/regression/"
AnalysisPath="/H4/regression/yasaman_sample/" # update here
setwd(path)
SummaryTable=read.table("SummaryTable.txt",header = TRUE,sep = ",")
#data index
#"RH.6hes_FPKM"      "B.P7.6hes_FPKM"    "B2.P11.6hes_FPKM"  "B2.P85.6hes_FPKM"  "B2.P148.6hes_FPKM" "B4.P12.6hes_FPKM"  "B4.P89.6hes_FPKM"

extraxCellulars=c(6,4,16,18) # extracellular samples only for B2

Intra=c(28, 2,12,14)

#response Label
ExtraSurvival6he=c(0.06109,0.093720,0.348810,0.166758709)
Invasion6he=c(0.1279182, 0.1411897,0.38589,0.3641795)
InvasionIntra=c(0.1279182, 0.1411897,0.3839,0.3670744)
PlauqSizeIntra=c(0.003615418,0.004486573,0.011716,0.013017)
PlauqSize6he=c(0.003615418,0.004486573,0.011793,0.013057)

Combined=c(PlauqSize6he,PlauqSizeIntra)

### start from zero
ExtraSurvival6he=c(0.06109,0.093720,0.348810,0.166758709)
Invasion6he=c(0.089767, 0.1411897,0.38589,0.3641795)
InvasionIntra=c(0.089767, 0.1411897,0.3839,0.3670744)
PlauqSizeIntra=c(0.002989,0.004486573,0.011716,0.013017)
PlauqSize6he=c(0.002989,0.004486573,0.011793,0.013057)

Combined=c(PlauqSize6he,PlauqSizeIntra)



# subsetting data

ExtradData=SummaryTable[,extraxCellulars]
rownames(ExtradData)=SummaryTable[,1]
head(ExtradData)

IntraData=SummaryTable[,Intra]
rownames(IntraData)=SummaryTable[,1]
head(IntraData)

CombinedData=SummaryTable[,c(extraxCellulars,Intra)]
rownames(CombinedData)=SummaryTable[,1]
head(CombinedData)
############################################
### Bootstrap create bootstrap directory
############################################
## create directory called bootstrap under analysispath then run the below code 

setwd(AnalysisPath)
i=1
tmp={}

BootstrapwithCV(t(ExtradData),PlauqSize6he ,1,Output=paste0(Output="bootstrap/","PlauqeSize_combined_"))

BootstrapwithCV(t(ExtradData), Invasion6he,1,Output=paste0(Output="bootstrap/","Invasion_combined_"))

BootstrapwithCV(t(ExtradData),ExtraSurvival6he ,1,Output=paste0(Output="bootstrap/","extraSurvival_combined_"))

#############################################
# combined table
#############################################
setwd(paste0(AnalysisPath,"/bootstrap/"))
dirname=list.dirs(path = ".")
nonzerodf=data.frame(condition="test",geneName="gene",coef=0,outSample="this")
for( i in 2:length(dirname))
    {
          fileList=list.files(path=dirname[i],pattern = "geneCoef.txt",full.names = T)
            conditioni=strsplit( dirname[i],split  ="./")[[1]][2]

            for(j in 1:length(fileList)){
                    cat(fileList[j],"\n")
                        sample=strsplit(fileList[j],split = "Sample_")[[1]][2]
                        sample=strsplit(sample,split = "_FPKM")[[1]][1]
                        nonzero=read.csv(fileList[j], header = T)
                        con=rep(conditioni,dim(nonzero)[1])
                        TMP=cbind(con,nonzero)

                        sample=rep(sample,dim(nonzero)[1])
                        TMP=cbind(TMP,sample)

                        colnames(TMP)=c("condition","geneName","coef","outSample")
                        nonzerodf=rbind(nonzerodf,TMP)
                      }


          }

nonzerodf=nonzerodf[-1,]
write.csv(nonzerodf,file="combinedNonzero.csv")


#############################################
# create the count of gene table
###################################################



conditions=as.character(unique(nonzerodf$condition))
out={}
for(i in 1:length(conditions))
    {
          tmp=nonzerodf[nonzerodf$condition==conditions[i],]
            a=table(tmp$geneName)
            index=which(a>0)
            write.csv(a[index],file=paste0(conditions[i],".csv"))
          }


########################
### feature selection
#########################

invatsionGenes=read.csv("Invasion_combined_1.csv")
plaqueGenes=read.csv("PlauqeSize_combined_1.csv")

ExtraSurvivalGenes=read.csv("extraSurvival_combined_1.csv")

sort(invatsionGenes$Freq,decreasing = T)
sort(plaqueGenes$Freq,decreasing = T)
sort(ExtraSurvivalGenes$Freq,decreasing = T)

o1=which(invatsionGenes$Freq>=45)
o2=which(plaqueGenes$Freq>=48)
o3=which(ExtraSurvivalGenes$Freq>=46)

selectedGeneNames=as.character(invatsionGenes$Var1[o1])
InvasionIndex={}
for(i in selectedGeneNames){
      InvasionIndex=c(InvasionIndex,which(rownames(ExtradData)==i))
      }

selectedGeneNames=as.character(plaqueGenes$Var1[o2])
plaqueIndex={}
for(i in selectedGeneNames){
      plaqueIndex=c(plaqueIndex,which(rownames(ExtradData)==i))
      }

selectedGeneNames=as.character(ExtraSurvivalGenes$Var1[o3])
ExtraSurvivalIndex={}
for(i in selectedGeneNames){
      ExtraSurvivalIndex=c(ExtraSurvivalIndex,which(rownames(ExtradData)==i))
      }


path="/H4/regression/yasaman_sample/withoutCV/"
setwd(paste0(AnalysisPath,"/withoutCV/"))
tmp={}
alpha=seq(0,1,0.05)
RSquareMatrix=matrix(0,nrow=1,ncol=3)

for(i in alpha){
      tmp={}
        output=LassoRegression(t(ExtradData[plaqueIndex,]), PlauqSize6he,i,Output=paste0(path,"PlauqeSize_"))
        tmp=c(tmp,output$rsq)



        output=LassoRegression(t(ExtradData[InvasionIndex,]),Invasion6he,i,Output=paste0(path,"Invasion_"))
        tmp=c(tmp,output$rsq)

        output=LassoRegression(t(ExtradData[ExtraSurvivalIndex,]),ExtraSurvival6he,i,Output=paste0(path,"extraSuerviva_"))
        tmp=c(tmp,output$rsq)
        RSquareMatrix=rbind(RSquareMatrix,tmp)
      }


RSquareMatrix=RSquareMatrix[-1,]
rownames(RSquareMatrix)=as.character(alpha)
colnames(RSquareMatrix)=c("Plauqe Size ","Invasion ", "Extra Survival")

det=(RSquareMatrix)

dataL = melt(det, id="x")


ggplot(dataL, aes(x=Var1, y=value, col=Var2))+
      geom <- line()+
      labs(title="R^2 Goodness of Fit across different alpha for B2 Samples without CV",x="Alpha",y="R^2")


labs(title="R^2 Goodness of Fit across different alpha for B2 Samples witht CV",x="Alpha",y="R^2")


dirname=list.dirs(path = ".")
nonzerodf=data.frame(condition="test",geneName="gene",coef=0,outSample="this")
for( i in 2:length(dirname))
    {
          fileList=list.files(path=dirname[i],pattern = "geneCoef.txt",full.names = T)
            conditioni=strsplit( dirname[i],split  ="./")[[1]][2]

            for(j in 1:length(fileList)){
                    cat(fileList[j],"\n")
                        sample=strsplit(fileList[j],split = "Sample_")[[1]][2]
                        sample=strsplit(sample,split = "_FPKM")[[1]][1]
                        nonzero=read.csv(fileList[j], header = T)
                        con=rep(conditioni,dim(nonzero)[1])
                        TMP=cbind(con,nonzero)

                        sample=rep(sample,dim(nonzero)[1])
                        TMP=cbind(TMP,sample)

                        colnames(TMP)=c("condition","geneName","coef","outSample")
                        nonzerodf=rbind(nonzerodf,TMP)
                      }


          }

nonzerodf=nonzerodf[-1,]
write.csv(nonzerodf,file="combinedNonzero.csv")

# create the count of gene table
conditions=as.character(unique(nonzerodf$condition))
out={}
for(i in 1:length(conditions))
    {
          tmp=nonzerodf[nonzerodf$condition==conditions[i],]
            a=table(tmp$geneName)
            index=which(a>0)
            write.csv(a[index],file=paste0(conditions[i],".csv"))
          }


