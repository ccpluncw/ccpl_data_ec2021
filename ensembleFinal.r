library(dplyr)
library(tidyr)
library("emmeans")

fitLL <- function(data, y, x, fileName) {
  df.tmp <- data.frame(x = data[,x], y = data[,y])

  dat.nls <- NULL
  dat.nls <- tryCatch(
      {
				nls(y ~ a + (1-a)/(1+exp(b-c*x)), data = df.tmp, start=c(a=0.2, b=4, c=1.5), nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE), algorithm = "port")

      },
      error=function(cond) {
          # Choose a return value in case of error
          return(NULL)
      },
      finally={
      }
  )

  if(is.null(dat.nls)) {
		dat.nls <- tryCatch(
      {
 				nls(y ~ a + (1-a)/(1+exp(b-c*x)), data = df.tmp, start=c(a=0.3, b=40, c=20), nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE), algorithm = "port")
	
      },
      error=function(cond) {
          # Choose a return value in case of error
          return(NULL)
      },
      finally={
      }
  	)
	}

  if(is.null(dat.nls)) {
		dat.nls <- tryCatch(
      {
 				nls(y ~ a + (1-a)/(1+exp(b-c*x)), data = df.tmp, start=c(a=0.4, b=.2, c=2), nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE), algorithm = "port")
	
      },
      error=function(cond) {
          # Choose a return value in case of error
          return(NULL)
      },
      finally={
      }
  	)
	}

  if(is.null(dat.nls)) {
		dat.nls <- tryCatch(
      {
 				nls(y ~ a + (1-a)/(1+exp(b-c*x)), data = df.tmp, start=c(a=0.4, b=.5, c=2), nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE), algorithm = "port")
	
      },
      error=function(cond) {
          # Choose a return value in case of error
          return(NULL)
      },
      finally={
      }
  	)
	}

		iterations <- 1
	  if(is.null(dat.nls)) {
			while(is.null(dat.nls)  & iterations < 1000 ) {
				a <- runif(1, min = -5, max = 5)
				b <- runif(1, min = -50, max = 50)
				c <- runif(1, min = -5, max = 5)
				dat.nls <- tryCatch(
		      {
		 				nls(y ~ a + (1-a)/(1+exp(b-c*x)), data = df.tmp, start=c(a=0.3, b=b, c=c), nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE), algorithm = "port")
	
		      },
		      error=function(cond) {
		          # Choose a return value in case of error
		          return(NULL)
		      },
		      finally={
		      }
		  	)
				iterations <- iterations + 1
			}
		}



  if(!is.null(dat.nls)) {
		df.tmp$Fit <- fitted(dat.nls)
    dat.r2.1 <- 1 - dat.nls$m$deviance()/sum((df.tmp$y - mean(df.tmp$y)^2))
    BIC.1 <- as.numeric(BIC(dat.nls))
    AICnorm.1 <- as.numeric(AIC(dat.nls))/length(df.tmp$y)
		
    a.1 <- as.numeric(coef(dat.nls)["a"])
    b.1 <- as.numeric(coef(dat.nls)["b"])
    c.1 <- as.numeric(coef(dat.nls)["c"])


    outList <- list(data = df.tmp, fit = dat.nls, r2 = dat.r2.1, BIC = BIC.1, AICnorm = AICnorm.1, a = a.1, b = b.1, c=c.1)

    pdf(fileName)
      with(df.tmp, plot(y ~ x))
      with(df.tmp[order(df.tmp$x),], lines(Fit ~ x, col="black"))
    dev.off()
  } else {
		print(fileName)
    outList <- list(data = df.tmp, fit = NA, r2 = NA, BIC = NA, AICnorm = NA, a = NA, b = NA, c=NA)

    pdf(fileName)
      with(df.tmp, plot(y ~ x))
    dev.off()
  }

  return(outList)
}

getLLdata <- function(data, pNoCol, totalCol, groupCol) {
  df.tmp <- data.frame(PNo = data[,pNoCol], tots = data[,totalCol],groupCol = data[,groupCol] )
  df.out <- data.frame(df.tmp %>% mutate(wPno = PNo*tots, .keep="all") %>% group_by_at(groupCol) %>% summarize(pNo = sum(wPno)/sum(tots), total = sum(tots)))
  colnames(df.out)[colnames(df.out) == "groupCol"] = groupCol
  return(df.out)
}

proportionPooledSE <- function(p1, p2, n1, n2) {
  # Calculate the pooled standard error
  pooled_se <- sqrt((p1*(1-p1)/n1) + (p2*(1-p2)/n2))
  return(pooled_se)
}

proportionSE <- function(p1, n1) {
  # Calculate the pooled standard error
  pSE <- sqrt(p1*(1-p1)/n1)
  return(pSE)
}

zTOp <- function(z, tail = "two") {
  # Check if tail is valid
  if (!tail %in% c("two", "one")) {
    stop("Invalid tail")
  }

  # Calculate the p-value
  if (tail == "two") {
    p_value <- 2 * pnorm(-abs(z))
  } else {
    p_value <- pnorm(z)
  }

  return(p_value)
}

getPredictedPoints <- function(data, pNoCol, distfrOldCol, distfrMeanCol, totalCol) {
		pNoN1A <- data[data[[distfrOldCol]] == 1 & data[[distfrMeanCol]] == 2,pNoCol]
		pNoN1B <- data[data[[distfrOldCol]] == -1 & data[[distfrMeanCol]] == -2,pNoCol]

		totN1A <- data[data[[distfrOldCol]] == 1 & data[[distfrMeanCol]] == 2,totalCol]
		totN1B <- data[data[[distfrOldCol]] == -1 & data[[distfrMeanCol]] == -2,totalCol]

		pNoN1 <- sum(pNoN1A*totN1A,pNoN1B*totN1B)/sum(totN1A,totN1B)
	
    #### This is predicted pNo assuming that a "different" is generated for both N1A and N1B
		pNoMeanSimilarity <- pNoN1A * pNoN1B
    pNoMeanActual <- data[data[[distfrMeanCol]] == 0,pNoCol]
    pNoOldActual <-  (data[data[[distfrOldCol]] == 0 & data[[distfrMeanCol]] == -1,pNoCol]+data[data[[distfrOldCol]] == 0 & data[[distfrMeanCol]] == 1,pNoCol])/2

    df.tmp <- data.frame(pNoMeanSimilarityPredicted = pNoMeanSimilarity, pNoMeanActual = pNoMeanActual, pNoOldActual =pNoOldActual, pNoN1A = pNoN1A, pNoN1B = pNoN1B, pNoN1 = pNoN1)

    return(df.tmp)
}

chanceThreshold <- 0.6
filename1 <- "mergedECdataRaw.txt"

df.raw.in <- read.table(filename1, header=T, sep="\t")
#### remove practice
df.raw.in <- df.raw.in[df.raw.in$Block > 1, ]
#### calculate absolute distance from mean)
df.raw.in$absDFM <- abs(df.raw.in$distfrMean)
df.raw.in$absDFO <- abs(df.raw.in$distfrOld)
### code NMO for relevant dimension
df.raw.in$relNMA <- ifelse(df.raw.in$absDFM == 0, "M", ifelse(df.raw.in$absDFO == 0, "O", "N"))


########################## Summarize Data ##########################

######### Here, we summarize 
df.raw.1 <- data.frame(df.raw.in %>% group_by(exp, sn, distfrMean, distfrOld, IrrDimLabel, distRelTto, IrrelDim, absDFM, absDFO, relNMA, Target, response) %>% summarize(n = length(Correct) ))
df.raw.2 <- df.raw.1 %>% pivot_wider(names_from = response, values_from = n)
df.raw.2$Y <- ifelse(is.na(df.raw.2$Y), 0, df.raw.2$Y)
df.raw.2$N <- ifelse(is.na(df.raw.2$N), 0, df.raw.2$N)
df.raw.2$tots <- df.raw.2$N + df.raw.2$Y
df.raw.2$PNo <- df.raw.2$N/df.raw.2$tots
df.raw <- df.raw.2

############# Summarize for irrelevant N1 only
#df.raw.n1.1 <- data.frame(df.raw.in[abs(df.raw.in$irrDFO) <2 , ] %>% group_by(exp, sn, distfrMean, distfrOld, IrrDimLabel2, distRelTto, IrrelDim, Target, response) %>% summarize(n = length(Correct) ))
df.raw.n1.1 <- data.frame(df.raw.in[abs(df.raw.in$irrDFO) <2 & abs(df.raw.in$distfrOld) <2 , ] %>% group_by(exp, sn, distfrMean, distfrOld, IrrDimLabel2, distRelTto, IrrelDim, absDFM, absDFO, relNMA,Target, response) %>% summarize(n = length(Correct) ))
df.raw.n1.2 <- df.raw.n1.1 %>% pivot_wider(names_from = response, values_from = n)
df.raw.n1.2$Y <- ifelse(is.na(df.raw.n1.2$Y), 0, df.raw.n1.2$Y)
df.raw.n1.2$N <- ifelse(is.na(df.raw.n1.2$N), 0, df.raw.n1.2$N)
df.raw.n1.2$tots <- df.raw.n1.2$N + df.raw.n1.2$Y
df.raw.n1.2$PNo <- df.raw.n1.2$N/df.raw.n1.2$tots
df.raw.n1 <- df.raw.n1.2


mainDir <- getwd()
newRunDir <- "analysis"
ch.newDir(mainDir, newRunDir)
analysisDir <- getwd()


exps <- unique(df.raw$exp)
for(ex in exps){
	df.raw.tmp <- df.raw[df.raw$exp == ex, ]
	df.raw.tmp.n1 <- df.raw.n1[df.raw.n1$exp == ex, ]
	
	condName <- gsub('/','',ex)
	condName <- gsub('-data','',condName)
	
	expDir <- paste(condName, "output", sep="_")
	ch.newDir(analysisDir, expDir)
	
	subs <- unique(df.raw.tmp$sn)

	df.llModelFits <- NULL
	df.simPointPredict <- NULL
	df.MNO <- NULL
	df.MNO.N1 <- NULL
	rm.sn <- NULL
	sn.keep <- NULL
	for(i in subs) {
		df.dat <- df.raw.tmp[df.raw.tmp$sn == i, ]
		df.dat.n1 <- df.raw.tmp.n1[df.raw.tmp.n1$sn == i, ]
	
		### analyze Distance From Mean
		### fit LL
		df.tmp.DFM <- getLLdata(df.dat, "PNo", "tots", "absDFM")
		if(df.tmp.DFM[6,"pNo"] > chanceThreshold) {
			fileName <- paste(i,"DFM.pdf", sep="-")
			dat.out.DFM <- fitLL(df.tmp.DFM, "pNo", "absDFM", fileName)
			DFM.ll <- as.data.frame(dat.out.DFM[3:8])
			DFM.ll$sn <- i
			DFM.ll$distRelTto <- df.dat$distRelTto[1]
			DFM.ll$IrrelDim <- df.dat$IrrelDim[1]
			DFM.ll$Target <- df.dat$Target[1]
		
			DFM.ll <- DFM.ll %>% rename(
				"DFM_r2" = "r2",
				"DFM_BIC" = "BIC",
				"DFM_AICnorm" = "AICnorm",
				"DFM_a" = "a",
				"DFM_b" = "b",
				"DFM_c" = "c"			
				)

			### analyze Distance From Old
			### fit LL
			df.tmp.DFO <- getLLdata(df.dat, "PNo", "tots", "absDFO")
			fileName <- paste(i,"DFO.pdf", sep="-")
			dat.out.DFO <- fitLL(df.tmp.DFO, "pNo", "absDFO", fileName)
			DFO.ll <- as.data.frame(dat.out.DFO[3:8])
			
			DFO.ll <- DFO.ll %>% rename(
				"DFO_r2" = "r2",
				"DFO_BIC" = "BIC",
				"DFO_AICnorm" = "AICnorm",
				"DFO_a" = "a",
				"DFO_b" = "b",
				"DFO_c" = "c"			
				)	
		
			df.ll.tmp <- cbind(DFM.ll, DFO.ll)
			df.llModelFits <- ch.rbind(df.llModelFits, df.ll.tmp)

			#predicted probability "No" to mean is the probability that both probes 1 point away say no
			df.tmp.DFO.1 <- data.frame(df.dat %>% mutate(wPno = PNo*tots, .keep="all") %>% group_by(distfrMean, distfrOld) %>% summarize(pNo = sum(wPno)/sum(tots), total = sum(tots)))
			df.tmp <- getPredictedPoints (df.tmp.DFO.1, "pNo", "distfrOld", "distfrMean", "total")
			df.tmp$sn <- i
			df.tmp$distRelTto <- df.dat$distRelTto[1]
			df.tmp$IrrelDim <- df.dat$IrrelDim[1]
			df.tmp$Target <- df.dat$Target[1]
	
			df.simPointPredict <- ch.rbind(df.simPointPredict, df.tmp)

			############################ NMO analysis ALL ############################
			df.tmp.NMO <- data.frame(df.dat %>% mutate(wPno = PNo*tots, .keep="all") %>% group_by(relNMA, IrrDimLabel) %>% summarize(pNo = sum(wPno)/sum(tots), total = sum(tots)))
	
			df.tmp.NMO$sn <- i
			df.tmp.NMO$distRelTto <- df.dat$distRelTto[1]
			df.tmp.NMO$IrrelDim <- df.dat$IrrelDim[1]
			df.tmp.NMO$Target <- df.dat$Target[1]
	
			df.MNO <- ch.rbind(df.MNO, df.tmp.NMO)

			############################ NMO analysis irrelevant N1 only ############################
			df.tmp.NMO.n1 <- data.frame(df.dat.n1 %>% mutate(wPno = PNo*tots, .keep="all") %>% group_by(relNMA, IrrDimLabel2) %>% summarize(pNo = sum(wPno)/sum(tots), total = sum(tots)))
	
			df.tmp.NMO.n1$sn <- i
			df.tmp.NMO.n1$distRelTto <- df.tmp.NMO.n1$distRelTto[1]
			df.tmp.NMO.n1$IrrelDim <- df.tmp.NMO.n1$IrrelDim[1]
			df.tmp.NMO.n1$Target <- df.tmp.NMO.n1$Target[1]
	
			df.MNO.N1 <- ch.rbind(df.MNO.N1, df.tmp.NMO.n1)

			############################ fill in above ############################
			

			sn.keep.tmp.M <- setNames(data.frame(t(df.tmp.DFM[,"pNo"])), df.tmp.DFM[,"absDFM"])
			sn.keep.tmp.M$sn <- i
			sn.keep.tmp.M$distRelTto <- df.dat$distRelTto[1]
			sn.keep.tmp.M$IrrelDim <- df.dat$IrrelDim[1]
			sn.keep.tmp.M$Target <- df.dat$Target[1]
			sn.keep.tmp.M$distance <- "DFM"
			sn.keep.tmp.O <- setNames(data.frame(t(df.tmp.DFO[,"pNo"])), df.tmp.DFO[,"absDFO"])
			sn.keep.tmp.O$sn <- i
			sn.keep.tmp.O$distRelTto <- df.dat$distRelTto[1]
			sn.keep.tmp.O$IrrelDim <- df.dat$IrrelDim[1]
			sn.keep.tmp.O$Target <- df.dat$Target[1]
			sn.keep.tmp.O$distance <- "DFO"
			sn.keep <- bind_rows(sn.keep, sn.keep.tmp.M, sn.keep.tmp.O)
	
	
		} else {
			rm.tmp <- setNames(data.frame(t(df.tmp.DFM[,"pNo"])), df.tmp.DFM[,"absDFM"])
			rm.tmp$sn <- i
			rm.tmp$distRelTto <- df.dat$distRelTto[1]
			rm.tmp$IrrelDim <- df.dat$IrrelDim[1]
			rm.tmp$Target <- df.dat$Target[1]

			rm.sn <- ch.rbind(rm.sn, rm.tmp)
		}
	
	}
	
	df.simPointPredict.long <- df.simPointPredict %>%
		pivot_longer(c("pNoMeanActual","pNoN1A","pNoN1B"), names_to = "trialType", values_to = "pNo")

	data1.RI.aov.n1rel=aov(pNo~trialType+Error(sn/(trialType)), data=df.simPointPredict.long)
 	options(contrasts = c("contr.sum", "contr.poly"))
 	emm.RI.relNMA.n1rel <- emmeans(data1.RI.aov.n1rel, ~ trialType)
 	emm.RI.relNMA.out.n1rel <- pairs(emm.RI.relNMA.n1rel)	

	
	ttestActVsN1 <- with(df.simPointPredict, t.test(pNoMeanActual, pNoN1, paired = T))
	ttestPredVsAct <- with(df.simPointPredict, t.test(pNoMeanActual, pNoMeanSimilarityPredicted, paired = T))
	ttestBIC <- with(df.llModelFits, t.test(DFM_BIC, DFO_BIC, paired = T))
	ttestAICnorm <- with(df.llModelFits, t.test(DFM_AICnorm, DFO_AICnorm, paired = T))
	ttestR2 <- with(df.llModelFits, t.test(DFM_r2, DFO_r2, paired = T))
	
	df.simPointPredict.sum <- df.simPointPredict %>% summarize(meanActualPno = mean(pNoMeanActual, na.rm = T), sdActualPno = sd(pNoMeanActual, na.rm = T), meanPredictPno = mean(pNoMeanSimilarityPredicted, na.rm = T), sdPredictPno = sd(pNoMeanSimilarityPredicted, na.rm = T), meanN1Pno = mean(pNoN1, na.rm = T), sdN1Pno = sd(pNoN1, na.rm = T))
	
	df.llModel.sum <- df.llModelFits %>% summarize(meanDFM_BIC = mean(DFM_BIC, na.rm = T), sdDFM_BIC = sd(DFM_BIC, na.rm = T), meanDFO_BIC = mean(DFO_BIC, na.rm = T), sdDFO_BIC = sd(DFO_BIC, na.rm = T), meanDFM_AICnorm = mean(DFM_AICnorm, na.rm = T), sdDFM_AICnorm = sd(DFM_AICnorm, na.rm = T), meanDFO_AICnorm = mean(DFO_AICnorm, na.rm = T), sdDFO_AICnorm = sd(DFO_AICnorm, na.rm = T), meanDFM_R2 = mean(DFM_r2, na.rm = T), sdDFM_R2 = sd(DFM_r2, na.rm = T), meanDFO_R2 = mean(DFO_r2, na.rm = T), sdDFO_R2 = sd(DFO_r2, na.rm = T))


#	df.MNO.N1$IrrDimLabel2 <- ifelse(df.MNO.N1$IrrDimLabel2 == "N1A" | df.MNO.N1$IrrDimLabel2 == "N1B", "N1", df.MNO.N1$IrrDimLabel2)


	if(df.MNO$IrrDimLabel[1] == "Na") {
		data1.RI.aov=aov(pNo~relNMA+Error(sn/(relNMA)), data=df.MNO)
 	 	options(contrasts = c("contr.sum", "contr.poly"))
 	 	emm.RI.relNMA <- emmeans(data1.RI.aov, ~ relNMA)
 	 	emm.RI.relNMA.out <- pairs(emm.RI.relNMA)	

		data1.RI.aov.n1=aov(pNo~relNMA+Error(sn/(relNMA)), data=df.MNO.N1)
		options(contrasts = c("contr.sum", "contr.poly"))
 	 	emm.RI.relNMA.n1 <- emmeans(data1.RI.aov.n1, ~ relNMA)
 	 	emm.RI.relNMA.out.n1 <- pairs(emm.RI.relNMA.n1)		
	} else {
		### run the MNO on all data
		data1.RI.aov=aov(pNo~relNMA*IrrDimLabel+Error(sn/(relNMA*IrrDimLabel)), data=df.MNO)
	
		 options(contrasts = c("contr.sum", "contr.poly"))
		 emm.RI.relNMA <- emmeans(data1.RI.aov, ~ relNMA)
		 emm.RI.relNMA.out <- pairs(emm.RI.relNMA)
		 emm.RI.IrrDimLabel <- emmeans(data1.RI.aov, ~ IrrDimLabel)
		 emm.RI.IrrDimLabel.out <- pairs(emm.RI.IrrDimLabel)

		 data1.R.M.aov=aov(pNo~IrrDimLabel+Error(sn/(IrrDimLabel)), data=df.MNO[df.MNO$relNMA == "M", ])
		 emm.R.M.IrrDimLabel <- emmeans(data1.R.M.aov, ~ IrrDimLabel)
		 emm.R.M.IrrDimLabel.out <- pairs(emm.R.M.IrrDimLabel)
		 data1.R.N.aov=aov(pNo~IrrDimLabel+Error(sn/(IrrDimLabel)), data=df.MNO[df.MNO$relNMA == "N", ])
		 emm.R.N.IrrDimLabel <- emmeans(data1.R.N.aov, ~ IrrDimLabel)
		 emm.R.N.IrrDimLabel.out <- pairs(emm.R.N.IrrDimLabel)
		 data1.R.O.aov=aov(pNo~IrrDimLabel+Error(sn/(IrrDimLabel)), data=df.MNO[df.MNO$relNMA == "O", ])
		 emm.R.O.IrrDimLabel <- emmeans(data1.R.O.aov, ~ IrrDimLabel)
		 emm.R.O.IrrDimLabel.out <- pairs(emm.R.O.IrrDimLabel)

 		### run the MNO on N1
 		data1.RI.aov.n1=aov(pNo~relNMA*IrrDimLabel2+Error(sn/(relNMA*IrrDimLabel2)), data=df.MNO.N1)
	
 		 options(contrasts = c("contr.sum", "contr.poly"))
 		 emm.RI.relNMA.n1 <- emmeans(data1.RI.aov.n1, ~ relNMA)
 		 emm.RI.relNMA.out.n1 <- pairs(emm.RI.relNMA.n1)
 		 emm.RI.IrrDimLabel.n1 <- emmeans(data1.RI.aov.n1, ~ IrrDimLabel2)
 		 emm.RI.IrrDimLabel.out.n1 <- pairs(emm.RI.IrrDimLabel.n1)

 		 data1.R.M.aov.n1=aov(pNo~IrrDimLabel2+Error(sn/(IrrDimLabel2)), data=df.MNO.N1[df.MNO.N1$relNMA == "M", ])
 		 emm.R.M.IrrDimLabel.n1 <- emmeans(data1.R.M.aov.n1, ~ IrrDimLabel2)
 		 emm.R.M.IrrDimLabel.out.n1 <- pairs(emm.R.M.IrrDimLabel.n1)
 		 data1.R.N.aov.n1=aov(pNo~IrrDimLabel2+Error(sn/(IrrDimLabel2)), data=df.MNO.N1[df.MNO.N1$relNMA == "N", ])
 		 emm.R.N.IrrDimLabel.n1 <- emmeans(data1.R.N.aov.n1, ~ IrrDimLabel2)
 		 emm.R.N.IrrDimLabel.out.n1 <- pairs(emm.R.N.IrrDimLabel.n1)
 		 data1.R.O.aov.n1=aov(pNo~IrrDimLabel2+Error(sn/(IrrDimLabel2)), data=df.MNO.N1[df.MNO.N1$relNMA == "O", ])
 		 emm.R.O.IrrDimLabel.n1 <- emmeans(data1.R.O.aov.n1, ~ IrrDimLabel2)
 		 emm.R.O.IrrDimLabel.out.n1 <- pairs(emm.R.O.IrrDimLabel.n1)

	 }

	 ##### create Figures
	 fileName <- paste("Data summary Plots.pdf", sep="-")
	 pdf(fileName)
		 # myColors <- c("grey90" , "grey70", "grey90", "grey90", "grey90")
		 #  		 with(df.simPointPredict, boxplot(pNoMeanSimilarityPredicted, pNoMeanActual, pNoOldActual, pNoN1A, pNoN1B, names = c("Similarity", "Mean", "Old", "N1A", "N1B"), col = myColors, ylim = c(0,1), frame = F, ylab = "p(No)"))
		 myColors <- c("grey90" , "grey70", "grey90", "grey90")
 		 with(df.simPointPredict, boxplot(pNoMeanSimilarityPredicted, pNoMeanActual, pNoOldActual, pNoN1, names = c("Similarity", "Mean", "Old", "N1"), col = myColors, ylim = c(0,1), frame = F, ylab = "p(No)"))
		 with(df.llModelFits, boxplot(DFM_BIC, DFO_BIC, names = c("DFM", "DFO"), ylim = c(-80,20), frame = F, ylab = "BIC"))	 
		 with(df.llModelFits, boxplot(DFM_AICnorm, DFO_AICnorm, names = c("DFM", "DFO"), ylim = c(-20,20), frame = F, ylab = "AIC norm"))	 
		 with(df.llModelFits, boxplot(DFM_r2, DFO_r2, names = c("DFM", "DFO"), ylim = c(.9, 1), frame = F, ylab = "r2"))
		 boxplot(sn.keep[sn.keep$distance == "DFM",c(1:6)], ylim=c(0,1), frame=F, ylab = "p(No)", main = "DFM")
		 boxplot(sn.keep[sn.keep$distance == "DFO",c(1:5)], ylim=c(0,1), frame=F, ylab = "p(No)", main = "DFO")	 

		 if(df.MNO.N1$IrrDimLabel2[1] != "Na") {
			df.MNO$IrrDimLabel <- factor(df.MNO$IrrDimLabel , levels=c("M", "O", "N"))
		 	with(df.MNO[df.MNO$relNMA == "M", ], boxplot(pNo ~ IrrDimLabel, names = c("Mean", "Old", "N"), ylim = c(0,1), frame = F, ylab = "p(No)", main = "All N, Relevant = Mean"))
		 	with(df.MNO[df.MNO$relNMA == "O", ], boxplot(pNo ~ IrrDimLabel, names = c("Mean", "Old", "N"), ylim = c(0,1), frame = F, ylab = "p(No)", main = "All N, Relevant = Old"))
		 	with(df.MNO[df.MNO$relNMA == "N", ], boxplot(pNo ~ IrrDimLabel, names = c("Mean", "Old", "N"), ylim = c(0,1), frame = F, ylab = "p(No)", main = "All N, Relevant = New"))
			df.MNO.N1$IrrDimLabel2 <- factor(df.MNO.N1$IrrDimLabel2 , levels=c("M", "O", "N1"))
		 	with(df.MNO.N1[df.MNO.N1$relNMA == "M", ], boxplot(pNo ~ IrrDimLabel2, names = c("Mean", "Old", "N1"), ylim = c(0,1), frame = F, ylab = "p(No)", main = "N1 Only, Relevant = Mean"))
		 	with(df.MNO.N1[df.MNO.N1$relNMA == "O", ], boxplot(pNo ~ IrrDimLabel2, names = c("Mean", "Old", "N1"), ylim = c(0,1), frame = F, ylab = "p(No)", main = "N1  Only, Relevant = Old"))
		 	with(df.MNO.N1[df.MNO.N1$relNMA == "N", ], boxplot(pNo ~ IrrDimLabel2, names = c("Mean", "Old", "N1"), ylim = c(0,1), frame = F, ylab = "p(No)", main = "N1  Only, Relevant = N1"))
		 }

	 dev.off()

	 filename <- paste(condName, "out.txt", sep="_")
	 sink(filename)
		 cat("\n Number of Subjects Removed for Chance Performance: ", nrow(rm.sn), " \n")
		 cat("\n removed subjects: \n")
		 print(rm.sn)
	 	
		 cat("\n\n ************************ ", condName, " ************************ \n\n")
		 cat("\n\t ******************** Similarity Hypothesis: Predicted vs Actual  ******************** \n\n")
		 writeLines(paste0("\t", capture.output(print(df.simPointPredict, row.names = F))))

		 cat("\n\t ********** Summary of Actual pNO vs Predicted pNo ********** \n\n")
		 writeLines(paste0("\t", capture.output(print(df.simPointPredict.sum, row.names = F))))
	
		 cat("\n\n\t T-test for pNO difference between Predicted vs Actual: \n")
		 writeLines(paste0("\t\t", capture.output(print(ttestPredVsAct, row.names = F))))

		 cat("\n\n\t T-test for pNO difference between N1 vs Actual: \n")
		 writeLines(paste0("\t\t", capture.output(print(ttestActVsN1, row.names = F))))

		 cat("\n\n\t Repeated Measures ANOVA for pNO difference between N1, N2, Actual: \n")
		 writeLines(paste0("\t\t", capture.output(print(summary(data1.RI.aov.n1rel), row.names = F))))
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.relNMA.n1rel, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.relNMA.out.n1rel, row.names = F))))
	
		 cat("\n\t ********** LL Model Fits: DFM vs DFO  ********** \n\n")
		 writeLines(paste0("\t", capture.output(print(df.llModelFits, row.names = F))))

		 cat("\n\t ********** Summary of LL Model Fits: DFM vs DFO  ********** \n\n")
		 writeLines(paste0("\t", capture.output(print(df.llModel.sum, row.names = F))))
	
		 cat("\n\n\t T-test for BIC difference between Predicted vs Actual: \n")
		 writeLines(paste0("\t\t", capture.output(print(ttestBIC, row.names = F))))
		 cat("\n\n\t T-test for AIC norm difference between Predicted vs Actual: \n")
		 writeLines(paste0("\t\t", capture.output(print(ttestAICnorm, row.names = F))))
		 cat("\n\n\t T-test for R2 difference between Predicted vs Actual: \n")
		 writeLines(paste0("\t\t", capture.output(print(ttestR2, row.names = F))))
	 
		 cat("\n\n ******************** M N O analysis on All Conditions: Repeated Measures Anova  ******************** \n\n")
		 writeLines(paste0("\t\t", capture.output(print(summary(data1.RI.aov)))))
		 cat("\n\n\t ********** Main Effects  ********** \n")
		 cat("\n\t\t ***** Target M N O  ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.relNMA, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.relNMA.out, row.names = F))))
	sink (NULL)
	
	if(df.MNO$IrrDimLabel[1] != "Na") {
		sink(filename, append =T)
		
		 cat("\n\n\t\t ***** Distractor M N O  ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.IrrDimLabel, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.IrrDimLabel.out, row.names = F))))
	
		 cat("\n\n\t ********** Interaction  ********** \n")
		 cat("\n\t\t ***** Distractor M N O  When Target = M ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.M.IrrDimLabel, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.M.IrrDimLabel.out, row.names = F))))
		 cat("\n\n\t\t ***** Distractor M N O  When Target = N ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.N.IrrDimLabel, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.N.IrrDimLabel.out, row.names = F))))
		 cat("\n\n\t\t ***** Distractor M N O  When Target = 0 ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print((emm.R.O.IrrDimLabel), row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.O.IrrDimLabel.out, row.names = F))))
	 	 
	 sink(NULL)
 }

	 ####
	sink(filename, append =T)
	 cat("\n\n ******************** M N O analysis on M O N1: Repeated Measures Anova  ******************** \n\n")
	 writeLines(paste0("\t\t", capture.output(print(summary(data1.RI.aov.n1)))))
	 cat("\n\n\t ********** Main Effects  ********** \n")
	 cat("\n\t\t ***** Target M N O  ***** \n")
	 cat("\n\t\t Means: \n\n")
	 writeLines(paste0("\t\t", capture.output(print(emm.RI.relNMA.n1, row.names = F))))
	 cat("\n\t\t post hoc: \n\n")
	 writeLines(paste0("\t\t", capture.output(print(emm.RI.relNMA.out.n1, row.names = F))))
	sink(NULL)

	if(df.MNO.N1$IrrDimLabel2[1] != "Na") {
		sink(filename, append =T)
		
		 cat("\n\n\t\t ***** Distractor M N1 O  ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.IrrDimLabel.n1, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.RI.IrrDimLabel.out.n1, row.names = F))))
	
		 cat("\n\n\t ********** Interaction  ********** \n")
		 cat("\n\t\t ***** Distractor M N1 O  When Target = M ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.M.IrrDimLabel.n1, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.M.IrrDimLabel.out.n1, row.names = F))))
		 cat("\n\n\t\t ***** Distractor M N1 O  When Target = N ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.N.IrrDimLabel.n1, row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.N.IrrDimLabel.out.n1, row.names = F))))
		 cat("\n\n\t\t ***** Distractor M N1 O  When Target = 0 ***** \n")
		 cat("\n\t\t Means: \n\n")
		 writeLines(paste0("\t\t", capture.output(print((emm.R.O.IrrDimLabel.n1), row.names = F))))
		 cat("\n\t\t post hoc: \n\n")
		 writeLines(paste0("\t\t", capture.output(print(emm.R.O.IrrDimLabel.out.n1, row.names = F))))
	 	 
	 sink(NULL)
	 
 }

	 setwd(analysisDir)
 }

 setwd(mainDir)

