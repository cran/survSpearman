## ----setup, include = TRUE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(survSpearman)
data(data123)

## ----ex10, echo = TRUE--------------------------------------------------------
data1 = data123[1:100,]
data1[1:7, ]

## ----exVis, echo = TRUE-------------------------------------------------------
visualBivarTimeToEvent(X = data1$X, Y = data1$Y, deltaX = data1$deltaX, deltaY = data1$deltaY,
                       xlim = c(0, 3), ylim = c(0, 3),
                       labelX = "X", labelY = "Y", segLength = 0.1, dotSize = 0.3,
                       scaleLegendGap = 1.1, legendCex = 0.6, labCex = 0.6, axisCex = 0.5)

## ----ex1a, echo = TRUE--------------------------------------------------------
res = survSpearman(X = data1$X, Y = data1$Y, deltaX = data1$deltaX, deltaY = data1$deltaY)
res[["Correlation"]]

## ----ex1b, echo = TRUE--------------------------------------------------------
dabrSurf = survDabrowska(data1$X, data1$Y, data1$deltaX, data1$deltaY)$DabrowskaEst
res = survSpearman(bivarSurf = dabrSurf, tauX = 2, tauY = 2)
res[["Correlation"]]

## ----ex11, echo = TRUE--------------------------------------------------------
dabrSurf = survDabrowska(data1$X, data1$Y, data1$deltaX, data1$deltaY)$DabrowskaEst
res = survSpearman(bivarSurf = dabrSurf, tauX = 2, tauY = 2)
res[["Restricted region set by user"]]
res[["Effective restricted region"]]

## ----ex12, echo = TRUE--------------------------------------------------------
res[["Correlation"]]

## ----ex13, echo = TRUE--------------------------------------------------------
res =  survSpearman(bivarSurf = dabrSurf, tauX = 1.5, tauY = 1.5)
res

## ----ex14, echo = TRUE--------------------------------------------------------
res =  survSpearman(bivarSurf = dabrSurf, tauX = 3, tauY = 3)
res

## ----ex20, echo = TRUE--------------------------------------------------------
data2 = data123[101:200,]
data2[data2$X == max(data2$X), c("X", "deltaX")]
data2[data2$Y == max(data2$Y), c("Y", "deltaY")]

## ----ex21, echo = TRUE--------------------------------------------------------
dabrSurf = survDabrowska(data2$X, data2$Y, data2$deltaX, data2$deltaY)$DabrowskaEst
res = survSpearman(bivarSurf = dabrSurf, tauX = Inf, tauY = Inf)
res

## ----ex22, echo = TRUE--------------------------------------------------------
data2[data2$X == max(data2$X), "deltaX"] = 0

## ----ex23, echo = TRUE--------------------------------------------------------
dabrSurf = survDabrowska(data2$X, data2$Y, data2$deltaX, data2$deltaY)$DabrowskaEst
res = survSpearman(bivarSurf = dabrSurf, tauX = Inf, tauY = Inf)
res

## ----ex30, echo = TRUE--------------------------------------------------------
data3 = data123[201:300,]
dabrSurf = survDabrowska(X = data3$X, Y = data3$Y, deltaX = data3$deltaX, deltaY = data3$deltaY)$DabrowskaEst
res = survSpearman(bivarSurf = dabrSurf, tauX = Inf, tauY = Inf)
res
cor(data3$X, data3$Y, method = "spearman")

## ----ex31, echo = TRUE--------------------------------------------------------
res =  survSpearman(bivarSurf = dabrSurf, tauX = qexp(0.75), tauY = qexp(0.75))
res

## ----ex40, echo = TRUE--------------------------------------------------------
data4 = data1[1:100,]
dabrSurf = survDabrowska(X = data4$X, Y = data4$Y, deltaX = data4$deltaX, deltaY = data4$deltaY)$DabrowskaEst
gridWidth = 0.1
survPMFPlot(bivarSurf = dabrSurf, gridWidthX = gridWidth, gridWidthY = gridWidth, scaleGapX = 1.5, scaleGapY = 1.5, XAxisLabel = "X", YAxisLabel = "Y", timeLabelScale = 1, axisLabelScale = 1, labelSkipX = 2, labelSkipY = 2, roundX = 2, roundY = 2, plotLabel = "Bivariate PMF")

## ----ex42, echo = TRUE--------------------------------------------------------
bivarSurfWithTails = rbind(cbind(dabrSurf, rep(0, nrow(dabrSurf))), rep(0, ncol(dabrSurf)+1))
lastRowVal = as.numeric(rownames(dabrSurf)[nrow(dabrSurf)]) + 2*gridWidth
lastColVal = as.numeric(colnames(dabrSurf)[ncol(dabrSurf)]) + 2*gridWidth
rownames(bivarSurfWithTails) = c(rownames(dabrSurf), as.character(lastRowVal))
colnames(bivarSurfWithTails) = c(colnames(dabrSurf), as.character(lastColVal))
survPMFPlot(bivarSurf = bivarSurfWithTails, gridWidthX = gridWidth, gridWidthY = gridWidth, scaleGapX = 1.5, scaleGapY = 1.5, XAxisLabel = "X", YAxisLabel = "Y", timeLabelScale = 1, axisLabelScale = 1, labelSkipX = 2, labelSkipY = 2, roundX = 2, roundY = 2, plotLabel = "Bivariate PMF")

## ----ex43, echo = TRUE--------------------------------------------------------
condSurf = survRestricted (bivarSurf = dabrSurf, tauX = 1.4, tauY = 1.4)$Sxy
survPMFPlot(bivarSurf = condSurf, gridWidthX = gridWidth, gridWidthY = gridWidth, scaleGapX = 1.5, scaleGapY = 1.5, XAxisLabel = "X", YAxisLabel = "Y", timeLabelScale = 1, axisLabelScale = 1, labelSkipX = 2, labelSkipY = 2, roundX = 2, roundY = 2, plotLabel = "Bivariate PMF")

## ----ex5a, echo = TRUE--------------------------------------------------------
### Histogram-like plot
plot(c(-0.5, 1.5), c(-0.5, 1), type = "n", xlab = "", ylab = "", cex.axis = 0.5)
vectValues = c(0.057, 0.336, 0.260, 0.362, 0.209)
plotVector(vectValues, width = 0.2, coordX = 0, coordY = 0, rotationRadians = 0,
   vectColor = "gray", bordColor = "white")

## ----ex5b1, echo = TRUE-------------------------------------------------------
width = c(0.10, 0.20, 0.10, 0.80, 0.12)
vectColor = c("orange", "green", "orchid", "blue", "goldenrod1")
plot(c(-0.5, 1.5), c(-0.5, 1.5), type = "n", xlab = "", ylab = "", cex.axis = 0.5)
plotVector(vectValues = vectValues, width, coordX = 0, coordY = 0,
   rotationRadians = pi/2, vectColor, bordColor = "white")
plotVector(vectValues = vectValues, width, coordX = 0, coordY = 0,
   rotationRadians = pi/4, vectColor, bordColor = "white")

## ----ex5b2, echo = TRUE-------------------------------------------------------
width = c(0.10, 0.20, 0.10, 0.80, 0.12)
vectColor = c("orange", "green", "orchid", "blue", "goldenrod1")
plot(c(-0.5, 1.5), c(-0.5, 1.5), type = "n", xlab = "", ylab = "", cex.axis = 0.5)
plotVector(vectValues = vectValues, width, coordX = 0, coordY = 0,
   rotationRadians = pi/2, vectColor, bordColor = "white")
plotVector(vectValues = -vectValues, width, coordX = 0, coordY = 0,
   rotationRadians = pi/2, vectColor, bordColor = "white")
plotVector(vectValues = -vectValues, width, coordX = 0, coordY = 0,
   rotationRadians = pi/4, vectColor, bordColor = "white")
plotVector(vectValues = -vectValues, width, coordX = 0, coordY = 0,
   rotationRadians = 0, vectColor, bordColor = "white")

## ----ex5c, echo = TRUE--------------------------------------------------------
vectValues = c(0.057, -0.336, 0.260, -0.222, 0.209)
plot(c(-1, 1), c(-0.5, 0.5), type = "n", xlab = "", ylab = "", cex.axis = 0.5)
vectColor = rep("goldenrod1", length(vectValues))
vectColor[vectValues<0] = "royalblue1"
bordColor = rep("red", length(vectValues))
bordColor[vectValues<0] = "darkblue"
plotVector(vectValues, width = 0.4, coordX = -1, coordY = 0, rotationRadians = 0,
    vectColor, bordColor = bordColor)

## ----ex6a, echo = TRUE--------------------------------------------------------
colorMatr = matrix(c("goldenrod1", "mediumpurple3", "palegreen3",
   "royalblue1", "orchid", "firebrick1"), nrow = 2, byrow = TRUE)
plot(c(1, 4), c(2, 4), type = "n", xlab = "", ylab = "", cex.axis = 0.5)
polygon(x = c(0, 5, 5, 0, 0), y = c(0, 0, 5, 5, 0), col = "black")
plotMatrix(colorMatr, borderCol = "black", coordX = 1, coordY = 2, widthX = 1, widthY = 1)

## ----ex6b, echo = TRUE--------------------------------------------------------
colorMatr = matrix(c("goldenrod1", "mediumpurple3", "palegreen3",
   "royalblue1", "orchid", "firebrick1"), nrow = 2, byrow = TRUE)
plot(c(1, 4.5), c(2, 4), type = "n", xlab = "", ylab = "", cex.axis = 0.5)
plotMatrix(colorMatr, borderCol = "white", coordX = 1, coordY = 2,
   widthX = c(1/4, 1, 1/6), widthY = c(1, 1/2))
reshapeColMatr = matrix(t(colorMatr), ncol = 1, nrow = nrow(colorMatr)*ncol(colorMatr),
   byrow = TRUE)
plotMatrix(reshapeColMatr, borderCol = "white", coordX = 2.8, coordY = 2,
   widthX = c(1/4), widthY = c(1/4))
text(x = rep(3.03, nrow(reshapeColMatr)),
   y = 2.12+c(0,cumsum(rep(1/4, nrow(reshapeColMatr)-1))),
   labels = reshapeColMatr[nrow(reshapeColMatr):1], pos = 4, cex = 0.5)

## ----ex7a, echo = TRUE--------------------------------------------------------
fillUpStr(c("A", "aaaa", NA, "bb", "4.5"), where="tail", fill = "-")

## ----ex7b, echo = TRUE--------------------------------------------------------
fillUpStr(c("A", "aaaa", NA, "bb", "4.5"), where="tail", howManyToFill = 10, fill = "-")
fillUpStr(c("!", "####", NA, "<<", "4.5"), where="head", howManyToFill = 10, fill = "*")

## ----ex8a, echo = TRUE--------------------------------------------------------
fillUpDecimals(c("2", "3.4", "A", NA))
fillUpDecimals(c("2", "3.4", "A.5", NA), howManyToFill = 5)
fillUpDecimals(c("2", "3.4", "A.xyz", NA), fill = "?")

## ----ex8a1, echo = TRUE-------------------------------------------------------
fillUpDecimals(c("2", "3.4", "A", NA), howManyToFill = -3)

## ----ex8b, echo = TRUE--------------------------------------------------------
fillUpDecimals(numVec = c("2", "3.4", "A.zx", NA, "#", "#a"), fill = "?")

