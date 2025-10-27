import argparse
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import cv2
import math
import json
import numpy as np
import matplotlib.pyplot as plt
from extPoint import find_extrema_points
from ccHull import create_convex_hull
from cmfPo import create_mask_from_polygon
from tranDueFuns import translate_points

def main():
    parser = argparse.ArgumentParser(description="Registration coordinate update.")
    parser.add_argument("-j", "--jsonFile", type=str, default='Adjusted_Output.json',
                        help="The JSON file with coordinates (default: Adjusted_Output.json)") 
    parser.add_argument("-cp", "--cellPointsFile", type=str, default='cellPoints.txt',
                        help="The file with cell points (default: cellPoints.txt)")
    parser.add_argument("-s", "--size", type=int, nargs=2, default=(23593, 22342),
                        help="The size of the original image in pixels, format: width height (default: 23593 22342)")
    parser.add_argument("-S", "--signal", action="store_true", help="Add Signal field if set") 
    
    args = parser.parse_args()  
    jsonFile = args.jsonFile    
    cellPointsFile = args.cellPointsFile
    size = args.size
    directory = os.path.dirname(jsonFile)

    if not directory:
        directory = os.getcwd()
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

    B = np.loadtxt(cellPointsFile, delimiter='\t', dtype=float)
    B = create_convex_hull(B)
    B = B.copy()

    maskB = create_mask_from_polygon(B, size)
    centerB = find_extrema_points(maskB)
    plt.imshow(maskB, cmap='gray')
    plt.axis('off')  
    cv2.imwrite(os.path.join(directory, 'cellPoints_MASK.png'), maskB)    

    with open(jsonFile, 'r') as file:
        data = json.load(file)

    shapes = data['shapes']
    
    combined_array = []
    if args.signal:
        sign_points = None
        for shape in shapes:
            if shape['label'] == 'Sign':
                sign_points = shape['points']
                break  

        if sign_points is not None:
            combined_array = sign_points
        else:
            print("No 'Sign' label found.")
    else:      
        for shape in shapes:
            if shape['label'] in ['Wai', 'Nei']:
                combined_array.extend(shape['points'])
    
    AN = create_convex_hull(np.array(combined_array))
    AN = AN.copy()
    maskA = create_mask_from_polygon(AN, size)

    plt.imshow(maskA, cmap='gray')
    plt.axis('off')
    updateFile = os.path.join(directory, 'markingPoints_MASK.png')
    cv2.imwrite(updateFile, maskA)

if __name__ == "__main__":
    main()

