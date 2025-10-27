rds = 'modified_data.rds'
json = 'SignalMatrix/AllRegistrationSchemes.json'
idjson = 'SignalMatrix/image_dimensions.json'  # 维度文件	
pixelssDNAFile = NULL
pixelHEFile = NULL
index = 2   # 选择哪一个索引对应的数据集	

library(SpLin)
RapidRegisterOUT(rds, json, idjson, pixelssDNAFile = pixelssDNAFile, pixelHEFile = pixelHEFile, index = index)

