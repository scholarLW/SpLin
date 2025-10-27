import argparse
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import csv
import cv2
import math
import numpy as np
import matplotlib.pyplot as plt
from extPoint import find_extrema_points
from ccHull import create_convex_hull
from cmfPo import create_mask_from_polygon
from tranDueFuns import translate_points

def main():
    parser = argparse.ArgumentParser(description="Registration coordinate update.")
    parser.add_argument("-mp", "--markingPointsFile", type=str, default='markingPoints.txt',
                        help="The file with marking points (default: markingPoints.txt)")
    parser.add_argument("-cp", "--cellPointsFile", type=str, default='cellPoints.txt',
                        help="The file with cell points (default: cellPoints.txt)")
    parser.add_argument("-s", "--size", type=int, nargs=2, default=(23593, 22342),
                        help="The size of the original image in pixels, format: width height (default: 23593 22342)")                     

    args = parser.parse_args()
    markingPointsFile = args.markingPointsFile
    cellPointsFile = args.cellPointsFile
    size = args.size

    directory = os.path.dirname(markingPointsFile)
    if not directory:
        directory = os.getcwd()
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

    AOLD = []
    with open(markingPointsFile, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:      
            AOLD.append(row)

    A = [[item['X'], item['Y']] for item in AOLD]
    A = np.array(A, dtype=float)
    A = create_convex_hull(A)
    A = A.copy()
    maskA = create_mask_from_polygon(A, size)
    centerA = find_extrema_points(maskA)
    plt.imshow(maskA, cmap='gray')
    plt.axis('off') 
    cv2.imwrite(os.path.join(directory, 'HE_pixel_MASK.png'), maskA)

    BOLD = []
    with open(cellPointsFile, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')  
        for row in reader:
            BOLD.append(row)

    B = [[item['X'], item['Y']] for item in BOLD]
    B = np.array(B, dtype=float)
    B = create_convex_hull(B)
    B = B.copy()
    maskB = create_mask_from_polygon(B, size)
    centerB = find_extrema_points(maskB)
    plt.imshow(maskB, cmap='gray')
    plt.axis('off')  
    cv2.imwrite(os.path.join(directory, 'ssDNA_pixel_MASK.png'), maskB)    
    
    adjusted_shapes = []
    for row in AOLD:
        x = float(row['X'])
        y = float(row['Y'])
        points = [x, y]
        translated_points = translate_points([points], centerB, centerA)
        row['X'], row['Y'] = translated_points[0]
        adjusted_shapes.append(row)

    A = [[item['X'], item['Y']] for item in adjusted_shapes]
    A = np.array(A, dtype=float)
    A = create_convex_hull(A)
    A = A.copy()
    maskA = create_mask_from_polygon(A, size)
    plt.imshow(maskA, cmap='gray')
    plt.axis('off')  
    cv2.imwrite(os.path.join(directory, 'HE_pixel_MASK.png'), maskA)

if __name__ == "__main__":
    main()

