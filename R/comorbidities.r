
#######################
##Deyo-Charlson index##
#######################

#Combining the above functions, deyo takes a 5 digit ICD-9-CM code and produces a list of 3 items:
#1. the total Charlson Score 
#2. a binary data frame of which comorbidites patients have - 1 for having it and 0 if not
#3. The point value data frame of comorbidites - 1-6 for having it, and 0 for not
deyo <- function(input.frame) {
	# Convert icd9 codes to numeric values and convert V434X to 44390 
	apply.icd9.deyo <- function(input.frame) {
		ICD9.5digit.deyo <- function(icd.code){
			if (is.na(icd.code)) {icd.code <- "00000"}
			icd9.3 <- substr(icd.code, 1, 3)
			icd9.4 <- substr(icd.code, 4, 4)
			icd9.5 <- substr(icd.code, 5, 5)
			if (icd9.4 == "X") {icd9.4 <- 0}
			if (icd9.5 == "X") {icd9.5 <- 0}
			icd9.result <- paste(icd9.3, icd9.4, icd9.5, sep = "")
			if (icd9.result == "V4340") {icd9.result <- 44390}
			return(as.numeric(icd9.result)/100)
		}

		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=n.cols)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				output.frame[i,j] <- ICD9.5digit.deyo(input.frame[i,j])
			}
		}
		return(output.frame)
	}

	# Deal with NA values

	apply.convert.na <- function(input.frame) {
		convert.na <- function(input.val) {
			if (is.na(input.val)) {input.val <- 0}
			output.val <- input.val
			return(output.val)
		}
		
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=n.cols)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				output.frame[i,j] <- convert.na(input.frame[i,j])
			}
		}
		return(output.frame)
	}

	# The following function develops a matrix with rows devoted to respondetns and each 
	# column a comorbidity.  The number in each column is the point value of having the comorbidity
	comorbidities.deyo <- function(input.frame) {
		#create lists of comorbidities
		mi <- c(seq(from=410, to=410.9, by=0.01), 412)
		chf <- c(seq(from=428, to=428.9, by=0.01))
		pvd <- c(443.9, 441, 441.9, 785.4) #v code v43.4 not included in this list
		cvd <- c(seq(from=430, to=438, by=0.01))
		dementia <- c(seq(from=290, to=290.9, by=0.01))
		copd <- c(seq(from=490, to=496, by=0.01), seq(from=500, to=505, by=0.01), 506.4)
		rheum <- c(710, 710.1, 710.4, seq(from =714, to=714.2, by=0.01), 714.81, 725)
		pud <- c(seq(from=531, to=534.9, by=0.01))
		mild.liver <- c(571.2, 571.5, 571.6, seq(from=571.4, to=571.49, by=0.01))
		dm <- c(seq(from=250,to=250.3,by=0.01), 250.7)
		dm.comp <- c(seq(from=250.4, to=250.6, by=0.01)) #2 point items start here
		plegia <- c(344.1, seq(from=342, to=342.9, by=0.01))
		renal <- c(seq(from=582, to=582.9, by=0.01), seq(from=583, to=583.7, by=0.01), 585, 586,seq(from=588, to=588.9, by=0.01))
		malignancy <- c(seq(from=140, to=172.9, by=0.01), seq(from=174, to=195.8, by=0.01), seq(from=200, to=208.9, by=0.01))
		severe.liver <- c(seq(from=572.2, to=572.8, by=0.01),seq(from=456, to=456.21, by=0.01)) # 3 point item
		mets <- c(seq(from=196, to=199.1, by=0.01)) # 6 point items
		hiv <- c(seq(from=42, to=44.93, by=0.01))
		
		deyo.list <- list(mi,chf,pvd,cvd,dementia,copd,rheum,pud,mild.liver,dm,dm.comp,plegia,renal,malignancy,severe.liver,mets,hiv)
		
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=17)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				for (k in 1:length(deyo.list)) {
					for (m in 1:length(deyo.list[[k]])) {
						if (input.frame[i, j] == deyo.list[[k]][m]) {
							output.frame[i,k] <- 1
						}
					}
				}
			}
		}
		
		output.frame <- as.data.frame(output.frame)
		colnames(output.frame) <- c("MI","CHF","PVD","CVD","DEMENTIA","COPD","RHEUM","PUD","MILD.LIVER","DM","DM.COMP","PLEGIA","RENAL","MALIGNANCY","SEVERE.LIVER", "METS", "HIV")
		return(output.frame)
	}

	# Convert the frame of point values to a frame of 0 for not having and 1 for having
	convert.to.points <- function(input.frame) {
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- input.frame
		output.frame[,11] <- output.frame[,11] *2
		output.frame[,12] <- output.frame[,12] *2
		output.frame[,13] <- output.frame[,13] *2
		output.frame[,14] <- output.frame[,14] *2
		output.frame[,15] <- output.frame[,15] *3
		output.frame[,16] <- output.frame[,17] *6
		output.frame[,16] <- output.frame[,17] *6
		return(output.frame)
	}

	#The following function sums the points in the comorbidites matrix produced above
	total.points <- function (input.frame) {
		n.rows <- length(input.frame[,1])
		output.vector <- matrix(0, nrow=n.rows, ncol=1)
		for (i in 1:n.rows) {
			output.vector[i] <- sum(input.frame[i,])
		}
		return(output.vector)
	}

		
	
	interim.frame.1 <- apply.icd9.deyo(input.frame)
	interim.frame.2 <- apply.convert.na(interim.frame.1)
	interim.frame.3 <- comorbidities.deyo(interim.frame.2)
	interim.frame.4 <- convert.to.points(interim.frame.3)
	POINTS <- total.points(interim.frame.4)
	deyo.data <- list(POINTS, interim.frame.3,interim.frame.4)
	names(deyo.data) <- c("CHARLSON.SCORE", "COMORBIDITIES", "COMORBIDITIES.POINTS")
	return(deyo.data)
}


#############################
##Original elixhauser index##
#############################

# elixhauser() takes a 5 digit ICD-9-CM code and produces a list of 2 items:
#1. the total count of elixhauser comorbidities 
#2. a binary data frame of which comorbidites patients have - 1 for having it and 0 if not
elixhauser <- function(input.frame) {
	# Convert icd9 codes to numeric values and convert v codes
	apply.icd9.elixhauser <- function(input.frame) { 
		ICD9.5digit.elixhauser <- function(icd.code){ 
			process.v.codes <- function(v.code) {
				icd9.2.5 <- as.numeric(substr(v.code, 2, 5))
				if (icd9.2.5 == 4500) {v.code <- 42610}
				if (icd9.2.5 == 5330) {v.code <- 42610}
				if (icd9.2.5 == 4220) {v.code <- 09320}
				if (icd9.2.5 == 4330) {v.code <- 09320}
				if (icd9.2.5 == 4340) {v.code <- 44000}
				if (icd9.2.5 == 4200) {v.code <- 40311}
				if (icd9.2.5 == 4510) {v.code <- 40311}
				if (icd9.2.5 == 5600) {v.code <- 40311}
				if (icd9.2.5 == 5680) {v.code <- 40311}
				if (icd9.2.5 == 4270) {v.code <- 07032}
				if (icd9.2.5 == 1271) {v.code <- 53170}
				if ((icd9.2.5 >= 1000) & (icd9.2.5 <= 1090)) {v.code <- 14000}
				if (icd9.2.5 == 1071) {v.code <- 20000}
				if (icd9.2.5 == 1072) {v.code <- 20000}
				if (icd9.2.5 == 1079) {v.code <- 20000}
				if (icd9.2.5 == 1130) {v.code <- 29110}
				
				return (v.code)
			}
			
			if (is.na(icd.code)) {icd.code <- "00000"}
			icd9.1 <- substr(icd.code, 1, 1)
			icd9.3 <- substr(icd.code, 1, 3)
			icd9.4 <- substr(icd.code, 4, 4)
			icd9.5 <- substr(icd.code, 5, 5)
			if (icd9.4 == "X") {icd9.4 <- 0}
			if (icd9.5 == "X") {icd9.5 <- 0}
			icd9.result <- paste(icd9.3, icd9.4, icd9.5, sep = "")
			if (icd9.1 == "V") {icd9.result <- process.v.codes(icd9.result)}
				
			return(as.numeric(icd9.result)/100)
		}
	
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=n.cols)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				output.frame[i,j] <- ICD9.5digit.elixhauser(input.frame[i,j])
			}
		}
		return(output.frame)
	}

	apply.convert.na <- function(input.frame) {
		convert.na <- function(input.val) {
			if (is.na(input.val)) {input.val <- 0}
			output.val <- input.val
			return(output.val)
		}
		
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=n.cols)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				output.frame[i,j] <- convert.na(input.frame[i,j])
			}
		}
		return(output.frame)
	}
	
	# The following function develops a matrix with rows devoted to respondents and each 
	# column a comorbidity.
	points.elixhauser.30 <- function(input.frame) {
		#create lists of comorbidities
		chf <- c(398.91,402.11,402.91,404.11,404.13,404.91,404.93,seq(from=428, to=428.9, by=0.01))
		arrhythmia <- c(426.1,426.11,426.13,seq(from=426.2, to=426.53, by=0.01),seq(from=426.6, to=426.89, by=0.01),427,427.2,427.31,427.6,427.9,785)
		valve <- c(seq(from=93.2, to=93.24, by=0.01),seq(from=394, to=397.1, by=0.01),seq(from=424, to=424.91, by=0.01),seq(from=746.3, to=746.6, by=0.01)) 
		pulm.circ <- c(seq(from=416, to=416.9, by=0.01), 417.9) 
		pvd <- c(seq(from=440, to=440.9, by=0.01),441.2,441.4,441.7,441.9,seq(from=443.1, to=443.9, by=0.01),447.1,557.1,557.9)
		htn <- c(401.1,401.9,402.1,402.9,404.1,404.9,405.11,405.19,405.91,405.99) 
		paralysis <- c(seq(from =342, to=342.12, by=0.01), seq(from=342.9, to=344.9, by=0.01)) 
		neuro.other <- c(331.9,332,333.4,333.5,seq(from=334, to=335.9, by=0.01),340,seq(from=341.1, to=341.9, by=0.01),seq(from=345, to=345.11, by=0.01),seq(from=345.4, to=345.51, by=0.01),seq(from=345.8, to=345.91, by=0.01),348.1,348.3,780.3,784.3)
		chronic.pulm <- c(seq(from=490, to=492.8, by=0.01), seq(from=493, to=493.91, by=0.01),494,seq(from=495, to=505, by=0.01),506.4)
		dm.uncomp <- c(seq(from=250,to=250.33,by=0.01))
		dm.comp <- c(seq(from=250.4, to=250.73, by=0.01),seq(from=250.9, to=250.93, by=0.01))
		hypothyroid <- c(seq(from=243, to=244.2, by=0.01),244.8,244.9)
		renal <- c(403.11,403.91,404.12,404.92,585,586)
		liver <- c(70.32,70.33,70.54,456,456.1,456.2,456.21,571,571.2,571.3,seq(from=571.4, to=571.49, by=0.01),571.5,571.6,571.8,571.9,572.3,572.8)
		pud <- c(531.7,531.9,532.7,532.9,533.7,533.9,534.7,534.9) 
		hiv <- c(seq(from=42, to=44.9, by=0.01)) 
		lymphoma <- c(seq(from=200,to=202.38, by=0.01),seq(from=202.5,to=203.01, by=0.01),seq(from=203.8,to=203.81, by=0.01),238.6,273.3)
		mets <- c(seq(from=196,to=199.1, by=0.01))
		solid.tumor <- c(seq(from=140,to=172.9, by=0.01),seq(from=174,to=175.9, by=0.01),seq(from=179,to=195.8, by=0.01))
		rheum <- c(701,seq(from=710,to=710.9, by=0.01),seq(from=714,to=714.9, by=0.01),seq(from=720,to=720.9, by=0.01),725)
		coag <- c(seq(from=286.0,to=286.9, by=0.01),287.1,seq(from=287.3,to=287.5, by=0.01))
		obesity <- c(278)
		wt.loss <- c(seq(from=260,to=263.9, by=0.01))
		lytes <- c(seq(from=276,to=276.9, by=0.01))
		anemia.loss <- c(280)
		anemia.def <- c(seq(from=280.1,to=281.9, by=0.01),285.9)
		etoh <- c(291.1,291.2,291.5,291.8,291.9,seq(from=303.9,to=303.93, by=0.01),seq(from=305,to=305.03, by=0.01))
		drugs <- c(292,seq(from=292.82,to=292.89, by=0.01),292.9,seq(from=304,to=304.93, by=0.01),seq(from=305.2,to=305.93, by=0.01))
		psychoses <- c(seq(from=295,to=298.9, by=0.01),seq(from=299.1,to=299.11, by=0.01))
		depression <- c(300.4,301.12,309,309.1,311)
		
		elixhauser.list <- list(chf,arrhythmia,valve,pulm.circ,pvd,htn,paralysis,neuro.other,chronic.pulm,dm.uncomp,dm.comp,hypothyroid,renal,liver,pud,hiv,lymphoma,mets,solid.tumor,rheum,coag,obesity,wt.loss,lytes,anemia.loss,anemia.def,etoh,drugs,psychoses,depression)
		
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=30)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				for (k in 1:length(elixhauser.list)){
					for (m in 1:length(elixhauser.list[[k]])) {
						if (input.frame[i, j] == elixhauser.list[[k]][m]) {
							output.frame[i,k] <- 1
						}
					}
				}
			}
		}
		
		#Apply the elixhauser hierarchy
		for (i in 1:length(output.frame[,1])){
			if (output.frame[i,11]==1) {output.frame[i,10] <- 0}
			if (output.frame[i,18]==1) {output.frame[i,19] <- 0}
		}
		
		output.frame <- as.data.frame(output.frame)
		colnames(output.frame) <- c("CHF","ARRHTHMIA","VALVE","PULM.CIRC","PVD","HTN","PARALYSIS","NEURO.OTHER","CHRONIC.PULM","DM.UNCOMP","DM.COMP","HYPOTHYROID","RENAL","LIVER","PUD","HIV","LYMPHOMA","METS","SOLID.TUMOR","RHEUM","COAG","OBESITY","WT.LOSS","LYTES","ANEMIA.LOSS","ANEMIA.DEF","ETOH","DRUGS","PSYCHOSES","DEPRESSION" )
		return(output.frame)
	}
	
	
#The following function sums the points in the comorbidites matrix produced above
	total.points <- function (input.frame) {
		n.rows <- length(input.frame[,1])
		output.vector <- matrix(0, nrow=n.rows, ncol=1)
		for (i in 1:n.rows) {
			output.vector[i] <- sum(input.frame[i,])
		}
		return(output.vector)
	}


	
	interim.frame.1 <- apply.icd9.elixhauser(input.frame)
	interim.frame.2 <- apply.convert.na(interim.frame.1)
	interim.frame.3 <- points.elixhauser.30(interim.frame.2)
	POINTS <- total.points(interim.frame.3)
	elixhauser.data <- list(POINTS, interim.frame.3)
	names(elixhauser.data) <- c("COMORBIDITY.CT", "COMORBIDITIES")
	return(elixhauser.data)
}


###################################################
##AHRQ comorbidites v3.6 index (newer elixhauser)##
###################################################

# ahrq_v3.6() takes a 5 digit ICD-9-CM code and produces a list of 2 items:
#1. the total count of elixhauser comorbidities 
#2. a binary data frame of which comorbidites patients have - 1 for having it and 0 if not
ahrq <- function(input.frame) {
	# Convert icd9 codes to numeric values and convert v codes
	apply.icd9.ahrq <- function(input.frame) { 
		ICD9.5digit.ahrq <- function(icd.code){ 
			process.v.codes <- function(v.code) {
				icd9.2.5 <- as.numeric(substr(v.code, 2, 5))
				#Valvular disease
				if (icd9.2.5 == 4220) {v.code <- 09320} 
				if (icd9.2.5 == 4330) {v.code <- 09320}
				#PVD
				if (icd9.2.5 == 4340) {v.code <- 44000} 
				#Renal Failure
				if (icd9.2.5 == 4200) {v.code <- 58530} 
				if (icd9.2.5 == 4510) {v.code <- 58530} 
				if ((icd9.2.5 >= 5600) & (icd9.2.5 <= 5632)) {v.code <- 58530}  
				if (icd9.2.5 == 5680) {v.code <- 58530} 
				if (icd9.2.5 == 4511) {v.code <- 58530}  
				if (icd9.2.5 == 4512) {v.code <- 58530}  
				#Liver Diseae
				if (icd9.2.5 == 4270) {v.code <- 07022}
				#Obsesity
				if ((icd9.2.5 >= 8530) & (icd9.2.5 <= 8539)) {v.code <- 02780}
				if ((icd9.2.5 >= 8541) & (icd9.2.5 <= 8545)) {v.code <- 02780}				
				if (icd9.2.5 == 8554) {v.code <- 02780}
				
				return (v.code)
			}
			
			if (is.na(icd.code)) {icd.code <- "00000"}
			icd9.1 <- substr(icd.code, 1, 1)
			icd9.3 <- substr(icd.code, 1, 3)
			icd9.4 <- substr(icd.code, 4, 4)
			icd9.5 <- substr(icd.code, 5, 5)
			if (icd9.4 == "X") {icd9.4 <- 0}
			if (icd9.5 == "X") {icd9.5 <- 0}
			icd9.result <- paste(icd9.3, icd9.4, icd9.5, sep = "")
			if (icd9.1 == "V") {icd9.result <- process.v.codes(icd9.result)}
				
			return(as.numeric(icd9.result)/100)
		}
	
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=n.cols)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				output.frame[i,j] <- ICD9.5digit.ahrq(input.frame[i,j])
			}
		}
		return(output.frame)
	}

	apply.convert.na <- function(input.frame) {
		convert.na <- function(input.val) {
			if (is.na(input.val)) {input.val <- 0}
			output.val <- input.val
			return(output.val)
		}
		
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=n.cols)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				output.frame[i,j] <- convert.na(input.frame[i,j])
			}
		}
		return(output.frame)
	}
	
	# The following function develops a matrix with rows devoted to respondents and each 
	# column a comorbidity.
	points.ahrq <- function(input.frame) {
		#create lists of comorbidities
		chf <- c(398.91,seq(from=428, to=428.9, by=0.01), 402.01,402.11,402.91, 404.01,404.11,404.91, 404.03,404.13,404.93) 
		valve <- c(seq(from=93.2, to=93.24, by=0.01),seq(from=394, to=397.1, by=0.01),397.9,seq(from=424, to=424.99, by=0.01),seq(from=746.3, to=746.6, by=0.01)) 
		pulm.circ <- c(seq(from =415.11, to=415.19, by=0.01),seq(from=416, to=416.9, by=0.01), 417.9) 
		pvd <- c(seq(from=440, to=440.9, by=0.01),seq(from=440.0, to=441.9, by=0.01),seq(from =442.0, to=442.9, by=0.01),seq(from =443.1, to=443.9, by=0.01),444.21,441.22,447.1,449.0,557.1,557.9)
		htn <- c(401.1,401.9,seq(from =642.00, to=642.04, by=0.01),401.0,437.2,seq(from =642.20, to=642.24, by=0.01),402.0,402.1,402.9,405.09,405.19,405.99, 402.01,402.11,402.91, 403.0,403.1,403.9,405.01,405.11,405.91,seq(from=642.10, to=642.14, by=0.01),403.01,403.11,403.91, 404.0,404.1,404.9, 404.01,404.11,404.91, 404.02,404.12,404.92, 404.03,404.13,404.93, seq(from =642.7, to=642.74, by=0.01),seq(from =642.9, to=642.94, by=0.01))  
		paralysis <- c(seq(from =342, to=344.9, by=0.01), seq(from=438.20, to=438.53, by=0.01),780.72) 
		neuro.other <- c(seq(from=330.0, to=331.9, by=0.01),332,333.4,333.5,333.7,333.71,333.72,333.79,333.85,333.94,seq(from=334.0, to=335.9, by=0.01),338.0,340.0,seq(from=341.1, to=341.9, by=0.01),seq(from=345.0, to=345.11, by=0.01),seq(from =345.2, to=345.3, by=0.01),seq(from=345.4, to=345.91, by=0.01),347.0,347.01,347.10,347.11,seq(from =649.4, to=649.44, by=0.01),786.7,seq(from =786.7, to=786.73, by=0.01),780.30,780.31,780.32,780.39,780.97,784.3)
		chronic.pulm <- c(seq(from=490, to=492.8, by=0.01), seq(from=493, to=493.92, by=0.01),seq(from =494.0, to=494.1, by=0.01),seq(from=495, to=505, by=0.01),506.4)
		dm.uncomp <- c(seq(from=250,to=250.33,by=0.01),seq(from=648.0, to=648.04, by=0.01),seq(from=249.0, to=249.31, by=0.01))
		dm.comp <- c(seq(from=250.4, to=250.93, by=0.01),775.1,seq(from=249.4, to=249.91, by=0.01))
		hypothyroid <- c(seq(from=243, to=244.2, by=0.01),244.8,244.9)
		renal <- c(585.3,585.4,585.5,585.6,585.9, 403.01,403.11,403.91,404.02,404.12,404.92, 404.03,404.13,404.93) 
		liver <- c(70.22,70.23,70.32,70.33,70.44,70.54,456,456.1,456.2,456.21,571,571.2,571.3,seq(from=571.4, to=571.49, by=0.01),571.5,571.6,571.8,571.9,572.3,572.8)
		pud <- c(531.41,531.51,531.61,531.7,531.71,531.91,532.41,532.51,532.61,532.7,532.71,532.91,533.41,533.51,533.61,533.7,533.71,533.91,534.41,534.51,534.61,534.7,534.71,534.91) 
		hiv <- c(seq(from=42, to=44.9, by=0.01)) 
		lymphoma <- c(seq(from=200,to=202.38, by=0.01),seq(from=202.5,to=203.01, by=0.01),238.6,273.3,seq(from=203.02,to=203.82, by=0.01))
		mets <- c(seq(from=196,to=199.1, by=0.01),seq(from=209.7,to=209.75, by=0.01),209.79,789.51)
		solid.tumor <- c(seq(from=140,to=172.9, by=0.01),seq(from=174,to=175.9, by=0.01),seq(from=179,to=195.8, by=0.01),seq(from=209.0,to=209.24, by=0.01),seq(from=209.25,to=209.3, by=0.01),seq(from=209.31,to=209.36, by=0.01),seq(from=258.01,to=258.03, by=0.01)  )
		rheum <- c(701,seq(from=710,to=710.9, by=0.01),seq(from=714,to=714.9, by=0.01),seq(from=720,to=720.9, by=0.01),725) 
		coag <- c(seq(from=286.0,to=286.9, by=0.01),287.1,seq(from=287.3,to=287.5, by=0.01),seq(from=649.3,to=649.34, by=0.01),289.84)
		obesity <- c(278.0,278.01,278.03,seq(from=649.1,to=649.14, by=0.01),793.91) 
		wt.loss <- c(seq(from=260,to=263.9, by=0.01),783.21,783.22)
		lytes <- c(seq(from=276,to=276.9, by=0.01))
		anemia.loss <- c(280,seq(from=648.2,to=648.24, by=0.01))
		anemia.def <- c(seq(from=280.1,to=281.9, by=0.01),seq(from=285.21,to=285.29, by=0.01),285.9) 
		etoh <- c(seq(from=291.0,to=291.3, by=0.01),291.5,291.8,291.81,291.82,291.89,291.9,seq(from=303.0,to=303.93, by=0.01),seq(from=305,to=305.03, by=0.01))
		drugs <- c(292,seq(from=292.82,to=292.89, by=0.01),292.9,seq(from=304,to=304.93, by=0.01),seq(from=305.2,to=305.93, by=0.01),seq(from=648.3,to=648.34, by=0.01))
		psychoses <- c(seq(from=295,to=298.9, by=0.01),299.1,299.11)
		depression <- c(300.4,301.12,309,309.1,311)
		
		ahrq.list <- list(chf,valve,pulm.circ,pvd,htn,paralysis,neuro.other,chronic.pulm,dm.uncomp,dm.comp,hypothyroid,renal,liver,pud,hiv,lymphoma,mets,solid.tumor,rheum,coag,obesity,wt.loss,lytes,anemia.loss,anemia.def,etoh,drugs,psychoses,depression)
		
		n.rows <- length(input.frame[,1])
		n.cols <- length(input.frame[1,])
		output.frame <- matrix(0, nrow=n.rows, ncol=29)
		for (i in 1:n.rows){
			for (j in 1:n.cols) {
				for (k in 1:length(ahrq.list)){
					for (m in 1:length(ahrq.list[[k]])) {
						if (input.frame[i, j] == ahrq.list[[k]][m]) {
							output.frame[i,k] <- 1
						}
					}
				}
			}
		}
		
		#Apply the elixhauser hierarchy
		for (i in 1:length(output.frame[,1])){
			if (output.frame[i,10]==1) {output.frame[i,9] <- 0}
			if (output.frame[i,17]==1) {output.frame[i,18] <- 0}
		}
		
		output.frame <- as.data.frame(output.frame)
		colnames(output.frame) <- c("CHF","VALVE","PULM.CIRC","PVD","HTN","PARALYSIS","NEURO.OTHER","CHRONIC.PULM","DM.UNCOMP","DM.COMP","HYPOTHYROID","RENAL","LIVER","PUD","HIV","LYMPHOMA","METS","SOLID.TUMOR","RHEUM","COAG","OBESITY","WT.LOSS","LYTES","ANEMIA.LOSS","ANEMIA.DEF","ETOH","DRUGS","PSYCHOSES","DEPRESSION" )
		return(output.frame)
	}

	#The following function sums the points in the comorbidites matrix produced above
	total.points <- function (input.frame) {
		n.rows <- length(input.frame[,1])
		output.vector <- matrix(0, nrow=n.rows, ncol=1)
		for (i in 1:n.rows) {
			output.vector[i] <- sum(input.frame[i,])
		}
		return(output.vector)
	}


	
	interim.frame.1 <- apply.icd9.ahrq(input.frame)
	interim.frame.2 <- apply.convert.na(interim.frame.1)
	interim.frame.3 <- points.ahrq(interim.frame.2)
	POINTS <- total.points(interim.frame.3)
	elixhauser.data <- list(POINTS, interim.frame.3)
	names(elixhauser.data) <- c("COMORBIDITY.CT", "COMORBIDITIES")
	return(elixhauser.data)
}
