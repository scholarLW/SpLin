rds = 'modified_data.rds'
image = 'SignalMatrix/signalMatrix.png'
json = 'SignalMatrix/signalMatrix.json' 
idjson = 'SignalMatrix/image_dimensions.json'

library(SpLin)
getPolygonPionts(rds, image, json, idjson, Signal = TRUE)
