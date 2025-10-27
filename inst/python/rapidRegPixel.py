import argparse
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import json
import cv2
import csv
import math
import copy
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from ccHull import create_convex_hull
from cmfPo import create_mask_from_polygon, create_mask_from_image
from extPoint import find_extrema_points
from calReCo import intersection as ITS
from calReCo import calculate_deviation_coefficients as caldc
from calReCo import calculate_registration_coefficient as calrc
from tranDueFuns import due_points, translate_points, rotate_points, rotate_point, get_convex_points, calculate_distance, calculate_rotation_angle, scale_polygon_points, scale_polygon_points_NEW


def reRapidRegPixel(pixelHEFile, data, size, index):
    directory = os.path.dirname(pixelHEFile)
    if not directory:
        directory = os.getcwd()

    directory = os.path.join(directory, 'Scheme', str(index))
    outputCoordinatesFile = os.path.join(directory, 'ssDNAHE_Adjusted_output_coordinates_and_colors.txt')

    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

    theta = data['THETA']
    ratios = data['RATIOS']
    transform = data['TRAN']
    centerB = data['CENT']
    if 'RATIOSUpdate' in data:
        RATIOSUpdate_values = data['RATIOSUpdate']
    else:
        RATIOSUpdate_values = None

    if 'THETAIN' in data:
        THETAIN_values = data['THETAIN']
    else:
        THETAIN_values = None        

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


def process_theta_dis(ratios, convex_pointsA, centerA, centerB, size, k, epsilon, maskB): 
    ANEW = []
    for x, y in convex_pointsA:
        points = [x, y]
        points = scale_polygon_points_NEW([points], ratios, centerA)
        ANEW.append(points[0])

    A = np.array(ANEW, dtype=float)
    maskA = create_mask_from_polygon(A, size)
    centerA = find_extrema_points(maskA)    
    
    ANEW1 = []
    for x, y in ANEW:
        points = [x, y]
        translated_points = translate_points([points], centerB, centerA)
        ANEW1.append(translated_points[0])

    A = np.array(ANEW1, dtype=float)
    maskA = create_mask_from_polygon(A, size)   
    d_array = caldc(maskB, maskA, k)
    rc, pards = calrc(d_array, epsilon=epsilon)
    return rc, pards, ratios, ANEW1


def process_dis(convex_pointsA, centerB, size, k, epsilon, ratios, maskB, thetaIN):
    if thetaIN != 0:
        convex_pointsB = []
        for x, y in convex_pointsA:
            points = [x, y]     
            points = rotate_point(points, math.radians(thetaIN), rotation_center = centerB)
            convex_pointsB.append(points)
    else:
        convex_pointsB = convex_pointsA

    ANEW = []
    for x, y in convex_pointsB:
        points = [x, y]
        points = scale_polygon_points_NEW([points], ratios, centerB)
        ANEW.append(points[0])
        
    A = np.array(ANEW, dtype=float)
    maskA = create_mask_from_polygon(A, size)
    d_array = caldc(maskB, maskA, k)
    rc, pards = calrc(d_array, epsilon=epsilon)
    return rc, pards, ANEW
   
def main():
    parser = argparse.ArgumentParser(description="Registration coordinate update.")       
    parser.add_argument("-p", "--pixelHEFile", type=str, default='HE/Adjusted_output_coordinates_and_colors.txt',
                        help="The file with pixel coordinates and color points of adjacent section HE staining (default: HE/Adjusted_output_coordinates_and_colors.txt)")
    parser.add_argument("-o", "--outputJsonFile", type=str, default='AllRegistrationSchemes.json',
                        help="The output JSON file with adjusted coordinates (default: AllRegistrationSchemes.json)")
    parser.add_argument("-mpi", "--markingPointsImage", type=str, default='markingPoints_MASK.png',
                        help="The marking points image file (default: markingPoints_MASK.png)")
    parser.add_argument("-mpj", "--markingPoints_MASKJSON", type=str, default='markingPoints_MASK.json',
                        help="The JSON file with marking points coordinates (default: markingPoints_MASK.json)")
    parser.add_argument("-cpi", "--cellPointsImage", type=str, default='cellPoints_MASK.png',
                        help="The cell points image file (default: cellPoints_MASK.png)")
    parser.add_argument("-cpj", "--cellPoints_MASKJSON", type=str, default='cellPoints_MASK.json',
                        help="The JSON file with cell points coordinates (default: cellPoints_MASK.json)")
    parser.add_argument("-s", "--size", type=int, nargs=2, default=(23593, 22342),
                        help="The size of the original image in pixels, format: width height (default: 23593 22342)")
    parser.add_argument("-t", "--transform", type=str, choices=['horizontal', 'vertical', None], default=None,
                        help="The transformation type, can be 'horizontal', 'vertical', or None (default: None)")   
    parser.add_argument("-th", "--thetaStep", type=float, default=0.5,
                        help="The step size for angle deviation iteration (default: 0.5)")
    parser.add_argument("-m", "--miuStep", type=float, default=0.01,
                        help="The step size for scaling deviation iteration (default: 0.01)")
    parser.add_argument("-n", "--Ngrid", type=int, default=3,
                        help="The number of grids composed of theta and miu, more value leads to more resource and time consumption (default: 3)")                         
    parser.add_argument("-k", "--kpixel", type=int, default=10,
                        help="Pixel interval for evaluating registration coefficient (default: 10)")
    parser.add_argument("-e", "--epsilon", type=float, default=0.001,
                        help="Error ratio caused by human lasso and annotation (default: 0.001)")
    parser.add_argument("--up", type=float, default=0,
                        help="Y-axis positive direction, positive value indicates stretching outside the concentric circle, negative value indicates scaling inside the concentric circle (default: 0)")
    parser.add_argument("--down", type=float, default=0,
                        help="Y-axis negative direction, positive value indicates stretching outside the concentric circle, negative value indicates scaling inside the concentric circle (default: 0)")
    parser.add_argument("--left", type=float, default=0,
                        help="X-axis negative direction, positive value indicates stretching outside the concentric circle, negative value indicates scaling inside the concentric circle (default: 0)")
    parser.add_argument("--right", type=float, default=0,
                        help="X-axis positive direction, positive value indicates stretching outside the concentric circle, negative value indicates scaling inside the concentric circle (default: 0)")
    parser.add_argument("--theta", type=float, default=0,
                        help="Rotation angle, positive value indicates counterclockwise rotation, negative value indicates clockwise rotation (default: 0)")
                 
    args = parser.parse_args()
    pixelHEFile = args.pixelHEFile
    outputJsonFile = args.outputJsonFile
    markingPointsImage = args.markingPointsImage
    markingPoints_MASKJSON = args.markingPoints_MASKJSON
    cellPointsImage = args.cellPointsImage
    cellPoints_MASKJSON = args.cellPoints_MASKJSON
    size = args.size
    transform = args.transform
    thetaStep = args.thetaStep
    miuStep = args.miuStep
    ngrid = args.Ngrid
    k = args.kpixel
    epsilon = args.epsilon
    up = args.up
    down = args.down
    left = args.left
    right = args.right
    thetaIN = args.theta

    directory = os.path.dirname(outputJsonFile)
    if not directory:
        directory = os.getcwd()
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

    with open(markingPoints_MASKJSON, 'r') as file:
        markingPoints_MASK = json.load(file)

    shapes = markingPoints_MASK['shapes']
    src_pts = []
    for shape in shapes:
        points = [[x, y] for x, y in shape['points']]
        src_pts.append(points)

    with open(cellPoints_MASKJSON, 'r') as file:
        cellPoints_MASK = json.load(file)

    shapes = cellPoints_MASK['shapes']
    dst_pts = []
    for shape in shapes:
        points = [[x, y] for x, y in shape['points']]
        dst_pts.append(points)

    A = np.array(src_pts[0])
    B = np.array(dst_pts[0])

    maskA = create_mask_from_image(markingPointsImage, transform=transform)
    maskB = create_mask_from_image(cellPointsImage, transform=None)
    
    centerA = find_extrema_points(maskA)
    centerB = find_extrema_points(maskB)
    
    if transform == 'vertical':
        A[:, 1] = 2 * centerA[1] - A[:, 1]
    elif transform == 'horizontal':
        A[:, 0] = 2 * centerA[0] - A[:, 0]    
    
    A = [due_points(point, centerA, centerB) for point in A]
    OP = [[x - centerB[0], y - centerB[1]] for x, y in A]
    OPNEW = [[x - centerB[0], y - centerB[1]] for x, y in B]
    thetas = [calculate_rotation_angle(op, opnew) for op, opnew in zip(OP, OPNEW)]
    theta0 = np.max(thetas) 
    dis_AA = [calculate_distance(point, centerA) for point in A]
    mean_dis_A = np.mean(dis_AA)    
    dis_BB = [calculate_distance(point, centerB) for point in B]
    mean_dis_B = np.mean(dis_BB)
    dis0 = mean_dis_B / mean_dis_A
    convex_pointsA = get_convex_points(maskA)
    ANEW = []
    for x, y in convex_pointsA:
        points = [x, y]     
        points = rotate_point(points, theta0, rotation_center = centerA)
        translated_points = translate_points([points], centerB, centerA)
        ANEW.append(translated_points[0])
    A = np.array(ANEW, dtype=float)
    maskAN = create_mask_from_polygon(A, size)
    centerA = find_extrema_points(maskAN)   
    data_to_store = []  
    dis_array = []    
    convex_pointsA = ANEW
    ratios1 = {
        'D': dis0, 
        'B': dis0, 
        'C': dis0,  
        'A': dis0  
    }
    
    dis_array.append(ratios1)
    results = Parallel(n_jobs=2)(delayed(process_theta_dis)(ratios, convex_pointsA, centerA, centerB, size, k, epsilon, maskB) for ratios in dis_array)

    RC, pards, ratios, upANEW = zip(*results)
    filtered_RC = [r for r in RC if r is not None]
    filtered_pards = [p for r, p in zip(RC, pards) if r is not None]
    filtered_ratios = [r_ for r, r_ in zip(RC, ratios) if r is not None]
    filtered_ANEW = [a for r, a in zip(RC, upANEW) if r is not None]

    for r, pas, ratio_dict, ANEWtmp in zip(filtered_RC, filtered_pards, filtered_ratios, filtered_ANEW):
        tmp = {
            'RC': r,
            'PARDS': pas,
            'THETA': theta0,
            'RATIOS': ratio_dict,
            'TRAN': transform,
            'MASKPOINT': deepcopy(ANEWtmp),
            'CENT': np.array(centerB, dtype=np.float32).tolist()
        }
        directoryD = os.path.join(directory, 'Scheme/1')
        FileD = os.path.join(directoryD, 'ssDNAHE_Adjusted_output_coordinates_and_colors.txt')
        reRapidRegPixel(pixelHEFile, tmp, size, 1)
        data_to_store.append(tmp)
        
        
    ratiosUpdata = {
        'D': 1 + left,    
        'B': 1 + right,    
        'C': 1 + down,    
        'A': 1 + up   
    }        

    data_to_storeNew = []
    if up or down or left or right or thetaIN:
        for item in data_to_store:
            item_copy = copy.deepcopy(item)
            rc, pards, ANEWNN = process_dis(item_copy['MASKPOINT'], item_copy['CENT'], size, k, epsilon, ratiosUpdata, maskB, thetaIN)
            tmpu = {
                'RC': rc,
                'PARDS': pards,
                'THETA': theta0,
                'THETAIN': thetaIN,
                'RATIOS': ratio_dict,
                'RATIOSUpdate': ratiosUpdata,
                'TRAN': transform,
                'MASKPOINT': ANEWNN,  
                'CENT': np.array(centerB, dtype=np.float32).tolist()               
            }
            reRapidRegPixel(pixelHEFile, tmpu, size, 2)
            data_to_storeNew.append(tmpu)

    data_to_storeF = data_to_store + data_to_storeNew        

    json_data = json.dumps(data_to_storeF, indent=4)  

    with open(outputJsonFile, 'w') as json_file:
        json_file.write(json_data)

    print(f'Data has been written to ({outputJsonFile})')


if __name__ == "__main__":
    main()

