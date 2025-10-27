json = 'SignalMatrix/Adjusted_signalMatrix.json'
cellPointsFile = 'SignalMatrix/cellPoints.txt' 
markingPointsImage = 'SignalMatrix/markingPoints_MASK.png'
markingPoints_MASKJSON = 'SignalMatrix/markingPoints_MASK.json'
cellPointsImage = 'SignalMatrix/cellPoints_MASK.png'
cellPoints_MASKJSON = 'SignalMatrix/cellPoints_MASK.json' 	
idjson = 'SignalMatrix/image_dimensions.json'
transform = NULL
kpixel = 10  
epsilon = 0.001 
interval = 0.05 
intervalAngle = 30 
up = 0.1
down = 0.1
left = -0.15
right = -0.15
theta = -10
spatialpointsize = 2

library(SpLin)
RapidRegister(json, cellPointsFile, markingPointsImage, markingPoints_MASKJSON, cellPointsImage, cellPoints_MASKJSON, idjson, transform = transform, kpixel = kpixel, epsilon = epsilon, interval = interval, intervalAngle = intervalAngle, up = up, down = down, left = left, right = right, theta = theta, spatialpointsize = spatialpointsize, Signal = TRUE)
