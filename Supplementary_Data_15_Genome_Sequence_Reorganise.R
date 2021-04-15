# Load libraries
library("Biostrings")
library("stringr")



################################################
###### Sequence Manipulation Functions    ######
################################################
updateRange <- function(destination, position, distance){
	if(str_count(destination, "_") == 1){
		return(destination)
	}else{
		destString <- str_split_fixed(destination, "_", 3)
		destNum <- as.numeric(destString[2:3])
		if(destNum[1] >= position){
			destNum[1] <- destNum[1] + distance
			if(destNum[2] != 999999999) destNum[2] <- destNum[2] + distance
			destString[2:3] <- str_pad(destNum, 9, pad="0")
		}
		return(paste(destString, collapse="_"))
	}
}

updateRangeCut <- function(destination, position, distance){
	if(str_count(destination, "_") == 1){
		return(destination)
	}else{
		destString <- str_split_fixed(destination, "_", 3)
		destNum <- as.numeric(destString[2:3])
		if(destNum[1] >= position){
			destNum[1] <- destNum[1] - distance
			if(destNum[2] != 999999999) destNum[2] <- destNum[2] - distance
			destString[2:3] <- str_pad(destNum, 9, pad="0")
		}
		return(paste(destString, collapse="_"))
	}
}


insertTransform <- function(from_to){
		distanceTb <- matrix(str_split_fixed(from_to, ":", 4), ncol=2, byrow=T)
		if(distanceTb[1,2] == ""){
			chrString <- str_extract(distanceTb[1, 1], "[A-C]\\d\\d")
			if(str_detect(distanceTb[1, 1], "top")) return(paste(c(chrString, str_pad(0, 9,pad = "0")), collapse="_")) else return(paste(c(chrString, str_pad(9, 9,pad = "9")), collapse="_"))
		}else{
			distanceTb1 <- str_split_fixed(distanceTb[,2], "_", 3)
			chrString <- unique(distanceTb1[, 1])
			distanceTb2 <- as.numeric(distanceTb1[, 2:3])
			distanceTb3 <- round(mean(distanceTb2[!distanceTb2%in%range(distanceTb2)]))
			distanceTb3 <- c(distanceTb3-1, distanceTb3)
			# cat(distanceTb3)
			return(paste(c(chrString, str_pad(distanceTb3, 9,pad = "0")), collapse="_"))
		}
}

moveRange <- function(from_to){
	# from_to <- paste(ABC_edits[i, c(5:8, 10:13)], collapse=":")
	temp <- str_split_fixed(from_to, ":", 8)
	from <- paste(temp[1:4], collapse=":")
	to <-  paste(temp[5:8], collapse=":")
	fromNum <-str_match(insertTransform(from), "\\_(\\d+)\\b")[,2]
	toNum <- str_match(insertTransform(to), "\\_(\\d+)\\b")[,2]
	chrString <- str_extract(from, "[A-C]\\d\\d")
	rangeReturn <- paste(c(chrString, fromNum, toNum), collapse="_")
	return(rangeReturn)
}

checkUpdate <- function(loc, mr, de){
	# if(order(c(mr, de, loc))[2] == 2) return(T) else return(F)
	allList <- c(mr, de, loc)
	rankList <- rank(allList)
	border <- rankList[1:2]
	remainRanks <- rankList[-c(1:2)]
	return(sapply(remainRanks, checkBetween, border))
}
checkBetween <- function(checkIn, border){
	return(sum(sign(checkIn - border))==0)
}



# Invert...
checkRange <- function(geneLoc){
	if(geneLoc =="") return(0) else{
		# geneLoc <- "A02_025619029_025619969"
		distanceTb1 <- str_split_fixed(geneLoc, "_", 3)
		# chr <- unique(distanceTb1[, 1])
		distanceTb2 <- as.numeric(distanceTb1[, 2:3])
		distanceTb3 <- abs(distanceTb2[2]- distanceTb2[1])
		return(distanceTb3)
	}
}


InvertFuncCorrect <- function(chr, from, to){
	# chr <-unlist(AllSeqs[ABC_inverts$V3[i]]); from <- paste(ABC_inverts[i, c(5:8)], collapse=":"); to <- paste(ABC_inverts[i, c(10:13)], collapse=":")
	invertFrom <- as.numeric(str_match(insertTransform(from), "\\_(\\d+)\\b")[,2])
	invertTo <- as.numeric(str_match(insertTransform(to), "\\_(\\d+)\\b")[,2])
	if(invertTo < invertFrom){temp <- invertFrom; invertFrom <- invertTo; invertTo <- temp}
	if(invertTo == 999999999) invertTo <- NA
	if(invertFrom == 0) invertFrom<- 1
	if(invertFrom == 1) part1 <- "" else part1 <- subseq(chr, 1, invertFrom-1)
	part2 <- subseq(chr, invertFrom, invertTo)
	if(is.na(invertTo)) part3 <- "" else part3 <- subseq(chr, invertTo+1)
	chrNew <- xscat(part1, reverseComplement(part2), part3)
	return(chrNew)
}
# end of invert

######################################################
###### End of Sequence Manipulation Functions   ######
######################################################





dateTag <- format(Sys.Date(),format="%Y-%m-%d")
ABC_edits <- fread("ABC_edits.csv", header=T)


logFile <- paste0("update_log_",dateTag,".txt")
cat("", file=logFile)
for(genome in LETTERS[1:3]){
	OriginalSeqs <- readDNAStringSet(paste0("sequences/",genome,"_afterInvert_",dateTag,".fa"))
	ABC_withUpdate <- ABC_edits[Operation != "Invert" & str_detect(order_loc, genome)]
	ABC_withUpdate$destination <- ABC_withUpdate[, insertTransform(paste(c(V15, V16, V17, V18), collapse=":")), by =1:nrow(ABC_withUpdate)]$V1
	ABC_withUpdate$MoveRange <- ABC_withUpdate[, moveRange(paste(c(V5, V6, V7, V8, V10, V11, V12, V13), sep=":")), by =1:nrow(ABC_withUpdate)]$V1
	# ABC_withUpdate[MoveRange != "", "MoveRange"] <- ABC_withUpdate[MoveRange != "", paste(V3, MoveRange, sep="_")]
	ABC_withUpdateDT <- ABC_withUpdate[, c("order_loc","Operation", "V3","destination", "MoveRange", "Orientation")]
	fwrite(ABC_withUpdateDT, file=paste0("ABC_withUpdateDT_",genome,".csv"))
	ABC_withUpdateDT$order <- 1:nrow(ABC_withUpdateDT)
	AllSeqs <- OriginalSeqs
	for(i in 1:nrow(ABC_withUpdate)){
		if(ABC_withUpdateDT$Operation[i] == "Insert"){
			destString <- str_split_fixed(ABC_withUpdateDT$destination[i], "_", 3)
			destNum <- as.numeric(destString[2:3])
			chr <- AllSeqs[destString[1]]
			insertPart <- unlist(AllSeqs[ABC_withUpdateDT$V3[i]])
			if(ABC_withUpdateDT$Orientation[i] == "rev") insertPart <- reverseComplement(insertPart)
			# insertSeq <- xscat(insertSeq, insertPart)
			if(is.na(destNum[2])){
				if(destNum[1] == 0) chrNew <- xscat(insertPart, chr) else chrNew <- xscat(chr, insertPart)
			}else{
				topPart <- subseq(chr, 1, destNum[1])
				bottomPart <- subseq(chr,destNum[2])
				chrNew <- xscat(topPart, insertPart, bottomPart)
			}
			AllSeqs[destString[1]] <- chrNew  #update chr

			# debug
			cat(i, "\t---", ABC_withUpdateDT$Operation[i], ABC_withUpdateDT$V3[i], "into", ABC_withUpdateDT$destination[i], "--- \n")
			cat("\t", destString[1], "Before", width(chr), "increase", length(insertPart), "Now", width(AllSeqs[destString[1]]), "\n")

			# update all the Move and desitnation entries
			if(any(str_detect(ABC_withUpdateDT$destination, destString[1])& ABC_withUpdateDT$order > i)) ABC_withUpdateDT[str_detect(destination, destString[1]) & order> i ,"destination"]<-
				sapply(unlist(ABC_withUpdateDT[str_detect(destination, destString[1])& order> i ,"destination"]), updateRange, destNum[1], length(insertPart))
			if(any(str_detect(ABC_withUpdateDT$MoveRange, destString[1])& ABC_withUpdateDT$order > i)) ABC_withUpdateDT[str_detect(MoveRange, destString[1])& order> i ,"MoveRange"]<- sapply(unlist(ABC_withUpdateDT[str_detect(MoveRange, destString[1])& order> i ,"MoveRange"]), updateRange, destNum[1], length(insertPart))
			# cat("--- After Insertion, Update related destination--- \n")
			# print(ABC_withUpdateDT[str_detect(destination, destString[1])&order > i, c("MoveRange", "destination")])
			# ABC_withUpdateDT[str_detect(destination, destString[1]) |str_detect(MoveRange, destString[1]), c("destination", "MoveRange")]
		}else{ #In case of Move
			# debug
			cat(i, "\t--- Move from", ABC_withUpdateDT$MoveRange[i], "into", ABC_withUpdateDT$destination[i], "--- \n")
			# Perform cutting
			moveString <- str_split_fixed(ABC_withUpdateDT$MoveRange[i], "_", 3)
			moveNum <- as.numeric(moveString[2:3])
			moveChr <- AllSeqs[moveString[1]]

			destString <- str_split_fixed(ABC_withUpdateDT$destination[i], "_", 3)
			destNum <- as.numeric(destString[2:3])
			destchr <- AllSeqs[destString[1]]

			if(moveString[1] == destString[1]){
				# cat("--- Move within the same chromosome, no need to update sequence--- \n")
				# Move within the same chromosome
				if(destNum[1] > moveNum[1]){
					# Move part is on the left of the desitination point
					if(moveNum[1] == 1){
						topPart <- ""
					}else{
						topPart <- subseq(moveChr, 1, moveNum[1] -1)
					}
					MovePart <- subseq(moveChr, moveNum[1], moveNum[2])
					if(destNum[1] == 999999999){
						betweenPart <- subseq(moveChr, moveNum[2] + 1)
						bottomPart <- ""
					}else{
						betweenPart <- subseq(moveChr, moveNum[2] + 1, destNum[1])
						bottomPart <- subseq(moveChr, destNum[2])
					}
					if(ABC_withUpdateDT$Orientation[i] == "rev") MovePart <- reverseComplement(MovePart)
					chrNew <- xscat(topPart, betweenPart, MovePart, bottomPart)
					# cat("--- Updating ",moveString[1],"from",width(AllSeqs[moveString[1]]), "to",width(chrNew),"--- \n")
					AllSeqs[moveString[1]] <- chrNew
				}else{
					if(destNum[1] == 0){
						topPart <- ""
						destNum[2] <- 1
					}else{
						topPart <- subseq(moveChr, 1, destNum[1])
					}
					betweenPart <- subseq(moveChr, destNum[2], moveNum[1] -1)
					if(moveNum[2] ==  999999999){
						MovePart <- subseq(moveChr, moveNum[1])
						bottomPart <- ""
					}else{
						MovePart <- subseq(moveChr, moveNum[1], moveNum[2])
						bottomPart <- subseq(moveChr, moveNum[2] + 1)
					}
					if(ABC_withUpdateDT$Orientation[i] == "rev") MovePart <- reverseComplement(MovePart)
					chrNew <- xscat(topPart, MovePart, betweenPart, bottomPart)
					# cat("--- Updating ",moveString[1],"from",width(AllSeqs[moveString[1]]), "to",width(chrNew),"--- \n")
					AllSeqs[moveString[1]] <- chrNew
				}
				# update all the Move and desitnation entries
				mr <- ABC_withUpdateDT$MoveRange[i]
				de <- ABC_withUpdateDT$destination[i]
				checkUpdateNeededMove <- checkUpdate(ABC_withUpdateDT$MoveRange,mr,de)& ABC_withUpdateDT$order > i
				if(any(checkUpdateNeededMove)){
					if(de > mr){
						ABC_withUpdateDT[checkUpdateNeededMove ,"MoveRange"] <- sapply(ABC_withUpdateDT$MoveRange[checkUpdateNeededMove], updateRangeCut, moveNum[1], width(MovePart))
					}else{
						ABC_withUpdateDT[checkUpdateNeededMove ,"MoveRange"] <- sapply(ABC_withUpdateDT$MoveRange[checkUpdateNeededMove], updateRange, moveNum[1], width(MovePart))
					}
				} 
				checkUpdateNeededDes <- checkUpdate(ABC_withUpdateDT$destination,mr,de)& ABC_withUpdateDT$order > i
				if(any(checkUpdateNeededDes)){
					if(de > mr){
						ABC_withUpdateDT[checkUpdateNeededDes ,"destination"] <- sapply(ABC_withUpdateDT$destination[checkUpdateNeededDes], updateRangeCut, moveNum[1], width(MovePart))
					}else{
						ABC_withUpdateDT[checkUpdateNeededDes ,"destination"] <- sapply(ABC_withUpdateDT$destination[checkUpdateNeededDes], updateRange, moveNum[1], width(MovePart))
					}
				}
			}else{
				# cat("--- Move between chromosomes--- \n")
				# cat(i, "\t--- Move from", moveString[1], "into", destString[1], "--- \n")
				# Moving between chromosomes
				# cut Moving part off
				# Move part is on the left of the desitination point
				if(moveNum[1] == 1){
					topPart <- ""
				}else{
					topPart <- subseq(moveChr, 1, moveNum[1] -1)
				}
				# MovePart <- subseq(moveChr, moveNum[1], moveNum[2])
				if(moveNum[2] ==  999999999){
					MovePart <- subseq(moveChr, moveNum[1])
					bottomPart <- ""
				}else{
					MovePart <- subseq(moveChr, moveNum[1], moveNum[2])
					bottomPart <- subseq(moveChr, moveNum[2] + 1)
				}
				cat("\t", moveString[1],"from",width(AllSeqs[moveString[1]]))
				AllSeqs[moveString[1]] <- xscat(topPart,bottomPart) #update Moving chr
				cat("to",width(AllSeqs[moveString[1]]),"--- \n")
				# update destination chr
				if(destNum[1] == 0){
					topPart <- ""
					bottomPart <- destchr
				}else{
					if(destNum[1] ==  999999999){
						topPart <- destchr
						bottomPart <- ""
					}else{
						topPart <- subseq(destchr, 1, destNum[1])
						bottomPart <- subseq(destchr, destNum[2])
					}
				}
				if(ABC_withUpdateDT$Orientation[i] == "rev") MovePart <- reverseComplement(MovePart)

				cat("\t",destString[1],"from",width(AllSeqs[destString[1]]))
				AllSeqs[destString[1]] <- xscat(topPart,MovePart,bottomPart) #update Moving chr
				cat("to",width(AllSeqs[destString[1]]),"--- \n")

				# update all the Move and desitnation entries
				if(any(str_detect(ABC_withUpdateDT$destination, destString[1])& ABC_withUpdateDT$order > i)) ABC_withUpdateDT[str_detect(destination, destString[1]) & order> i ,"destination"]<-
					sapply(unlist(ABC_withUpdateDT[str_detect(destination, destString[1])& order> i ,"destination"]), updateRange, destNum[1], width(MovePart))
				if(any(str_detect(ABC_withUpdateDT$MoveRange, destString[1])& ABC_withUpdateDT$order > i)) ABC_withUpdateDT[str_detect(MoveRange, destString[1])& order> i ,"MoveRange"]<-
					sapply(unlist(ABC_withUpdateDT[str_detect(MoveRange, destString[1])& order> i ,"MoveRange"]), updateRange, destNum[1], width(MovePart))

				# update all the Move and desitnation entries
				if(any(str_detect(ABC_withUpdateDT$destination, moveString[1])& ABC_withUpdateDT$order > i)) ABC_withUpdateDT[str_detect(destination, moveString[1]) & order> i ,"destination"]<-
					sapply(unlist(ABC_withUpdateDT[str_detect(destination, moveString[1])& order> i ,"destination"]), updateRangeCut, moveNum[1], width(MovePart))
				if(any(str_detect(ABC_withUpdateDT$MoveRange, moveString[1])& ABC_withUpdateDT$order > i)) ABC_withUpdateDT[str_detect(MoveRange, moveString[1])& order> i ,"MoveRange"]<-
					sapply(unlist(ABC_withUpdateDT[str_detect(MoveRange, moveString[1])& order> i ,"MoveRange"]), updateRangeCut, moveNum[1], width(MovePart))
			}
		}
		print(ABC_withUpdateDT[order > i, c("MoveRange", "destination")])
	}
	NewSeqs <- AllSeqs[!names(AllSeqs) %in% ABC_withUpdateDT[Operation=="Insert"]$V3]
	writeXStringSet(NewSeqs, file=paste0("sequences/",genome,"_afterUpdate",dateTag,".fa"))
	writeXStringSet(NewSeqs[grep("[A-C]\\d\\d", names(NewSeqs))], file=paste0("sequences/",genome,"_afterUpdate",dateTag,"_chr.fa"))
}
