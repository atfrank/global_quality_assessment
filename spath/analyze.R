setwd("~/Documents/GitSoftware/global_quality_assessment/spath/")
source("sensivity_analysis.R")

f <-combine_features_shifts_all()
#imp_l <- get_importance_by_group(f, method="larmord")
#imp_r <- get_importance_by_group(f, method="ramsey")
imp_c <- get_importance_by_group(f, method="consensus", verbose = TRUE)
#save(imp_l, imp_r, imp_c, file = "sensitivity.RData")
stop()

#imp <- imp_l
#mat <- subset(imp, nucleus=="C1'")$`%IncMSE`
#mat <- cbind(mat, subset(imp, nucleus=="C2'")$`%IncMSE`)
#mat <- cbind(mat, subset(imp, nucleus=="C3'")$`%IncMSE`)
#mat <- cbind(mat, subset(imp, nucleus=="C4'")$`%IncMSE`)
#mat <- cbind(mat, subset(imp, nucleus=="C5'")$`%IncMSE`)
#mat <- cbind(mat, subset(imp, nucleus=="C2")$`%IncMSE`)
#mat <- cbind(mat, subset(imp, nucleus=="C5")$`%IncMSE`)
#mat <- cbind(mat, subset(imp, nucleus=="C6")$`%IncMSE`)
#mat <- cbind(mat, subset(imp, nucleus=="C8")$`%IncMSE`)
#names <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")
#l_names <- c("hydrogen bonding","magnetic anisotropy","bond polarization","close contact","α","β","γ","δ", "ε","ζ","χ","ν0","ν1","ν2","ν3","ν4","stacking","ring current")
#barplot(0.1*sqrt(mat), beside = T, names.arg = names, legend.text = l_names, col=rainbow(n=18), ylab="sensivity (ppm)", xlab="nucleus", cex.axis = 1.0, cex.names = 1.0, ylim=c(0, 3.0),args.legend = list(x = "topleft", bty = "n", inset=c(0.3, -0.4)))
#


# plots the distribution off a given parameter

cols <- rainbow(5)
data <- combine_features_shifts("1R2P", "2LPS")
tmp <- subset(data$larmord, 
plot(density(tmp$rc), xlim =c(-1,2), col=cols[1], lwd="2")

data <- combine_features_shifts("2FRL", "2M22")
tmp <- data$larmord
lines(density(tmp$rc), col = cols[2], lwd="2")

data <- combine_features_shifts("2H2X", "2M21")
tmp <- data$larmord
lines(density(tmp$rc), col = cols[3], lwd="2")

data <- combine_features_shifts("2KFC", "2L1V")
tmp <- data$larmord
lines(density(tmp$rc), col = cols[4], lwd="2")


# get data that includes the resid number and model get data
# for a given nucleus determine the fluctuations in the ring current features over the trajectories
# make boxplots over the various test cases

data <- combine_features_shifts_all(include_these = c("resid", "model","id"))
tmp <- data$larmord


test2 <- ddply(.data = tmp, .variables = c("id","resid", "nucleus"), .fun = get_fluctuations, feature = "beta", circular_data = TRUE)
boxplot(sd~id, data = test2)

test2 <- ddply(.data = tmp, .variables = c("id","resid", "nucleus"), .fun = get_difference, feature = "rc", circular_data = FALSE)
boxplot(df~id, data = test2)

test2 <- ddply(.data = tmp, .variables = c("id","resid", "nucleus"), .fun = get_difference, feature = "beta", circular_data = TRUE)
boxplot(df~id, data = test2)
