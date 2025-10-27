rds = 'SignalMatrix/Update_Manual_modified_data.rds'
json = 'SignalMatrix/Adjusted_Output_Manual.json'
nlim = 0
drop = TRUE  
triangleMerge = TRUE 
triangle.probs = 0.2 

library(SpLin)
coordLinearEXP(rds, json, nlim = nlim, drop = drop, triangleMerge = triangleMerge, triangle.probs = triangle.probs)


