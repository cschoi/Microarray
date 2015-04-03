getwd()
setwd("~/Study/BioNetwork/GSE28521")

total_data <- read.table("norm_total_GSE28521.csv", sep =",", row.names = 1, header = T)

head(total_data)
dim(total_data)

total_group <- total_data[,c(4:35)]
total_group <- as.matrix(total_group)
head(total_group)


l_total <- c(rep(1, 16), rep(2, 16))
l_total
g <- factor(l_total)
g


##### 6. SAM analysis #####
library(samr)     					# SAM package loading

sg = sub("0","2",l_total)  				# grouping
sm= list(x = total_group, y= sg, logged2=TRUE)			# SAM input matrix
st = samr(sm,resp.type="Two class unpaired",nperm=100)	# sam test
dt = samr.compute.delta.table(st)			# delta table 작성 Computes tables of thresholds, cutpoints and corresponding False Discovery rates for SAM analysis
dt
d = 0.44         

pdf("GSE28521_samrPlot_total.pdf", height=10, width=15)
samr.plot(st,d)						#SAM plot
dev.off()

st = samr.compute.siggenes.table(st,d,sm,dt)		#result table 작성

names(st)									#element 확인 
head(st$genes.up)					#수행결과 확인

st$ngenes.up
st$ngenes.lo

genes_up <- st$genes.up
genes_down <- st$genes.lo

ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), function(y) y[slot])

genes_up[,2] <- ss(genes_up[,2],"g", 2)
head(genes_up)
IDrow <- as.numeric(genes_up[,2])
head(IDrow)

data_up <- total_data[IDrow,]
head(data_up)
head(genes_up)

Fold_qval <- genes_up[,7:8]
head(Fold_qval)
Gene_info <- data_up[,1:3]
head(Gene_info)

sam_up <- cbind(Gene_info, Fold_qval)
head(sam_up)

colnames(sam_up) <- c("ID", "Gene_title", "Gene_symbol", "Fold_change", "q_value")
head(sam_up)
str(sam_up)

sam_up$Gene_symbol <- as.character(paste(sam_up$Gene_symbol))
sam_up$q_value <- as.numeric(paste(sam_up$q_value))
sam_up$Fold_change <- as.numeric(paste(sam_up$Fold_change))

which(sam_up$q_value < 0.05)

write.table(sam_up,"sam_up.csv",sep=",", row.names=T)	#결과 파일로 쓰기

