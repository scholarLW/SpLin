json = 'SignalMatrix/Adjusted_signalMatrix.json'
cellPointsFile = 'SignalMatrix/cellPoints.txt'
idjson = 'SignalMatrix/image_dimensions.json' 

library(SpLin)
preAutoRegister(json, cellPointsFile, idjson, Signal = TRUE)
