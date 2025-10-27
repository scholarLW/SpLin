import argparse
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import json
import cv2
import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from ccHull import create_convex_hull
from cmfPo import create_mask_from_polygon, create_mask_from_image
from extPoint import find_extrema_points
from calReCo import intersection as ITS
from calReCo import calculate_deviation_coefficients as caldc
from calReCo import calculate_registration_coefficient as calrc
from tranDueFuns import due_points, translate_points, rotate_points, rotate_point, get_convex_points, \
    calculate_distance, calculate_rotation_angle, scale_polygon_points, scale_polygon_points_NEW


def main():
    parser = argparse.ArgumentParser(description="Registration coordinate update.")
    parser.add_argument("-HE", "--pixelHEFile", type=str, default=None,
                        help="The file with pixel coordinates and color points of HE staining (default: None")
    parser.add_argument("-ssDNA", "--pixelssDNAFile", type=str, default=None,
                        help="The file with pixel coordinates and color points of ssDNA staining (default: None")
    parser.add_argument("-j", "--jsonFile", type=str, default='AllRegistrationSchemes.json',
                        help="The JSON file with coordinates (default: AllRegistrationSchemes.json)")
    parser.add_argument("-s", "--size", type=int, nargs=2, default=(23593, 22342),
                        help="The size of the original image in pixels, format: width height (default: 23593 22342)")
    parser.add_argument("-i", "--index", type=int, default=1,
                        help="The index of the dataset to be selected (default: 1)")

    args = parser.parse_args()
    pixelHEFile = args.pixelHEFile
    pixelssDNAFile = args.pixelssDNAFile
    jsonFile = args.jsonFile
    index = args.index
    size = args.size

    directory = os.path.dirname(jsonFile)
    if not directory:
        directory = os.getcwd()
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    with open(jsonFile, 'r') as file:
        data = json.load(file)

    A = data[index - 1]['MASKPOINT']
    theta = data[index - 1]['THETA']
    ratios = data[index - 1]['RATIOS']
    transform = data[index - 1]['TRAN']
    centerB = data[index - 1]['CENT']
    jsonData = data[index - 1]['JSON']
    if 'RATIOSUpdate' in data[index - 1]:
        RATIOSUpdate_values = data[index - 1]['RATIOSUpdate']
    else:
        RATIOSUpdate_values = None

    if 'THETAIN' in data[index - 1]:
        THETAIN_values = data[index - 1]['THETAIN']
    else:
        THETAIN_values = None    
        
    jsonO = os.path.join(directory, 'Adjusted_Output_Manual.json')
    json_data = json.dumps(jsonData, indent=4) 
    with open(jsonO, 'w') as json_file:
        json_file.write(json_data)

    A = np.array(A, dtype=float)
    maskA = create_mask_from_polygon(A, size)
    plt.imshow(maskA, cmap='gray')
    plt.axis('off') 

    updateFile = os.path.join(directory, 'markingPoints_MASK_update.png')
    cv2.imwrite(updateFile, maskA)


    if pixelHEFile is not None:
        directory = os.path.dirname(pixelHEFile)
        if not directory:
            directory = os.getcwd()
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

        outputCoordinatesFile = os.path.join(directory, 'pixel_markingPoints_XY.txt')
        AOLDHE = []
        with open(pixelHEFile, newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t') 
            for row in reader:
                AOLDHE.append(row)

        AP = [[item['X'], item['Y']] for item in AOLDHE]
        AP = np.array(AP, dtype=float)
        AP = create_convex_hull(AP)
        AP = AP.copy()
        maskAP = create_mask_from_polygon(AP, size)
        centerAP = find_extrema_points(maskAP)

        plt.imshow(maskAP, cmap='gray')
        plt.axis('off')  
        cv2.imwrite(os.path.join(directory, 'HE_pixel_MASK_pixel.png'), maskAP)

        adjusted_shapes = []
        for row in AOLDHE:
            x = float(row['X'])
            y = float(row['Y'])
            points = [x, y]
            points = rotate_points(points, centerAP, transform=transform) 
            points = rotate_point(points, theta, rotation_center=centerAP)
            translated_points = translate_points([points], centerB, centerAP)
            row['X'], row['Y'] = translated_points[0]
            adjusted_shapes.append(row)

        AP = [[item['X'], item['Y']] for item in adjusted_shapes]
        AP = np.array(AP, dtype=float)
        APc = create_convex_hull(AP)
        APc = APc.copy()
        maskAP = create_mask_from_polygon(APc, size)
        centerAP = find_extrema_points(maskAP)

        adjusted_shapesA = []
        for row in adjusted_shapes:
            x = float(row['X'])
            y = float(row['Y'])
            points = [x, y]
            points = scale_polygon_points_NEW([points], ratios, centerAP)
            row['X'], row['Y'] = points[0]
            adjusted_shapesA.append(row)

        AP = [[item['X'], item['Y']] for item in adjusted_shapesA]
        AP = np.array(AP, dtype=float)
        AP = create_convex_hull(AP)
        AP = AP.copy()
        maskAP = create_mask_from_polygon(AP, size)
        centerAP = find_extrema_points(maskAP)

        adjusted_shapesB = []
        for row in adjusted_shapesA:
            x = float(row['X'])
            y = float(row['Y'])
            points = [x, y]
            translated_points = translate_points([points], centerB, centerAP)
            row['X'], row['Y'] = translated_points[0]
            adjusted_shapesB.append(row)

        AP = [[item['X'], item['Y']] for item in adjusted_shapesB]
        AP = np.array(AP, dtype=float)
        AP = create_convex_hull(AP)
        AP = AP.copy()
        maskAP = create_mask_from_polygon(AP, size)
        centerAP = find_extrema_points(maskAP)


        if THETAIN_values is not None and THETAIN_values != 0:
            adjusted_shapesC = []
            for row in adjusted_shapesB:
                x = float(row['X'])
                y = float(row['Y'])
                points = [x, y]
                points = rotate_point(points, math.radians(THETAIN_values), rotation_center = centerB)
                row['X'], row['Y'] = points
                adjusted_shapesC.append(row)
            AP = [[item['X'], item['Y']] for item in adjusted_shapesC]
            AP = np.array(AP, dtype=float)
            AP = create_convex_hull(AP)
            AP = AP.copy()
            maskAP = create_mask_from_polygon(AP, size)
        else:
            adjusted_shapesC = adjusted_shapesB

        if RATIOSUpdate_values is not None:
            adjusted_shapesD = []
            for row in adjusted_shapesC:
                x = float(row['X'])
                y = float(row['Y'])
                points = [x, y]
                points = scale_polygon_points_NEW([points], RATIOSUpdate_values, centerB)
                row['X'], row['Y'] = points[0]
                adjusted_shapesD.append(row)
            AP = [[item['X'], item['Y']] for item in adjusted_shapesD]
            AP = np.array(AP, dtype=float)
            AP = create_convex_hull(AP)
            AP = AP.copy()
            maskAP = create_mask_from_polygon(AP, size)
        else:
            adjusted_shapesD = adjusted_shapesC

        plt.imshow(maskAP, cmap='gray')
        plt.axis('off')  

        cv2.imwrite(os.path.join(directory, 'HE_pixel_MASK_pixel_Update.png'), maskAP)

        with open(outputCoordinatesFile, "w") as f:
            f.write("X\tY\tColor\n")
            for item in adjusted_shapesD:
                line = f"{item['X']:.6f}\t{item['Y']:.6f}\t{item['Color']}\n"
                f.write(line)



    if pixelssDNAFile is not None:
        directory = os.path.dirname(pixelssDNAFile)
        if not directory:
            directory = os.getcwd()
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

        outputCoordinatesFile = os.path.join(directory, 'pixel_markingPoints_XY.txt')

        AOLD = []
        with open(pixelssDNAFile, newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t') 
            for row in reader:
                AOLD.append(row)

        AP = [[item['X'], item['Y']] for item in AOLD]
        AP = np.array(AP, dtype=float)
        APc = create_convex_hull(AP)
        APc = APc.copy()
        maskAP = create_mask_from_polygon(APc, size)
        centerAP = find_extrema_points(maskAP)

        adjusted_shapes = []
        for row in AOLD:
            x = float(row['X'])
            y = float(row['Y'])
            points = [x, y]
            points = rotate_points(points, centerAP, transform=transform)  
            points = rotate_point(points, theta, rotation_center=centerAP)
            translated_points = translate_points([points], centerB, centerAP)
            row['X'], row['Y'] = translated_points[0]
            adjusted_shapes.append(row)

        AP = [[item['X'], item['Y']] for item in adjusted_shapes]
        AP = np.array(AP, dtype=float)
        APc = create_convex_hull(AP)
        APc = APc.copy()
        maskAP = create_mask_from_polygon(APc, size)
        centerAP = find_extrema_points(maskAP)

        adjusted_shapesA = []
        for row in adjusted_shapes:
            x = float(row['X'])
            y = float(row['Y'])
            points = [x, y]
            points = scale_polygon_points_NEW([points], ratios, centerAP)
            row['X'], row['Y'] = points[0]
            adjusted_shapesA.append(row)

        AP = [[item['X'], item['Y']] for item in adjusted_shapesA]
        AP = np.array(AP, dtype=float)
        AP = create_convex_hull(AP)
        AP = AP.copy()
        maskAP = create_mask_from_polygon(AP, size)
        centerAP = find_extrema_points(maskAP)

        adjusted_shapesB = []
        for row in adjusted_shapesA:
            x = float(row['X'])
            y = float(row['Y'])
            points = [x, y]
            translated_points = translate_points([points], centerB, centerAP)
            row['X'], row['Y'] = translated_points[0]
            adjusted_shapesB.append(row)

        AP = [[item['X'], item['Y']] for item in adjusted_shapesB]
        AP = np.array(AP, dtype=float)
        AP = create_convex_hull(AP)
        AP = AP.copy()
        maskAP = create_mask_from_polygon(AP, size)
        centerAP = find_extrema_points(maskAP)

        if THETAIN_values is not None and THETAIN_values != 0:
            # print(THETAIN_values)
            adjusted_shapesC = []
            for row in adjusted_shapesB:
                x = float(row['X'])
                y = float(row['Y'])
                points = [x, y]
                points = rotate_point(points, math.radians(THETAIN_values), rotation_center = centerB)
                row['X'], row['Y'] = points
                adjusted_shapesC.append(row)
            AP = [[item['X'], item['Y']] for item in adjusted_shapesC]
            AP = np.array(AP, dtype=float)
            AP = create_convex_hull(AP)
            AP = AP.copy()
            maskAP = create_mask_from_polygon(AP, size)
        else:
            adjusted_shapesC = adjusted_shapesB

        if RATIOSUpdate_values is not None:
            #print(RATIOSUpdate_values)
            adjusted_shapesD = []
            for row in adjusted_shapesC:
                x = float(row['X'])
                y = float(row['Y'])
                points = [x, y]
                points = scale_polygon_points_NEW([points], RATIOSUpdate_values, centerB)
                row['X'], row['Y'] = points[0]
                adjusted_shapesD.append(row)
            AP = [[item['X'], item['Y']] for item in adjusted_shapesD]
            AP = np.array(AP, dtype=float)
            AP = create_convex_hull(AP)
            AP = AP.copy()
            maskAP = create_mask_from_polygon(AP, size)
        else:
            adjusted_shapesD = adjusted_shapesC

        plt.imshow(maskAP, cmap='gray')
        plt.axis('off')  

        cv2.imwrite(os.path.join(directory, 'ssDNA_pixel_MASK_pixel_Update.png'), maskAP)

        with open(outputCoordinatesFile, "w") as f:
            f.write("X\tY\tColor\n")
            for item in adjusted_shapesD:
                line = f"{item['X']:.6f}\t{item['Y']:.6f}\t{item['Color']}\n"
                f.write(line)

    print("The coordinates have been updated and saved.")


if __name__ == "__main__":
    main()

