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

def process_theta_dis(ratios, shapes, centerA, centerB, size, k, epsilon, maskB, signal):
    adjusted_shapes1 = []   
    for shape in shapes:
        points = [[x, y] for x, y in shape['points']]
        points = scale_polygon_points_NEW(points, ratios, centerA)  
        shape['points'] = points
        adjusted_shapes1.append(shape)
    combined_array = []
    if signal:
        sign_points = None
        for shape in adjusted_shapes1:
            if shape['label'] == 'Sign':
                sign_points = shape['points']
                break 

        if sign_points is not None:
            combined_array = sign_points
        else:
            print("No 'Sign' label found.")
    else:      
        for shape in adjusted_shapes1:
            if shape['label'] in ['Wai', 'Nei']:
                combined_array.extend(shape['points'])
    
    AN1 = create_convex_hull(np.array(combined_array))
    maskAN1 = create_mask_from_polygon(AN1, size)
    centerAN1 = find_extrema_points(maskAN1)

    adjusted_shapes2 = []   
    for shape in adjusted_shapes1:
        points = [[x, y] for x, y in shape['points']]
        translated_points = translate_points(points, centerB, centerAN1)
        shape['points'] = translated_points
        adjusted_shapes2.append(shape)

    combined_array = []
    if signal:
        sign_points = None
        for shape in adjusted_shapes2:
            if shape['label'] == 'Sign':
                sign_points = shape['points']
                break 

        if sign_points is not None:
            combined_array = sign_points
        else:
            print("No 'Sign' label found.")
    else:     
        for shape in adjusted_shapes2:
            if shape['label'] in ['Wai', 'Nei']:
                combined_array.extend(shape['points'])    
                
    AN1 = create_convex_hull(np.array(combined_array))
    maskAN1 = create_mask_from_polygon(AN1, size)    
    d_array = caldc(maskB, maskAN1, k)
    rc, pards = calrc(d_array, epsilon=epsilon)
    return rc, pards, ratios, adjusted_shapes2, AN1

def process_dis(data, centerB, size, k, epsilon, ratios, maskB, thetaIN, signal):
    shapes = data['shapes']
    if thetaIN != 0:
        adjusted_shapesA = []        
        for shape in shapes:
            points = [rotate_point([x, y], math.radians(thetaIN), rotation_center = centerB) for x, y in shape['points']]
            shape['points'] = points
            adjusted_shapesA.append(shape)            
    else:
        adjusted_shapesA = shapes

    adjusted_shapes = []   
    for shape in adjusted_shapesA:
        points = [[x, y] for x, y in shape['points']]
        points = scale_polygon_points_NEW(points, ratios, centerB)
        shape['points'] = points
        adjusted_shapes.append(shape)
    combined_array = []
    if signal:
        sign_points = None
        for shape in adjusted_shapes:
            if shape['label'] == 'Sign':
                sign_points = shape['points']
                break 

        if sign_points is not None:
            combined_array = sign_points
        else:
            print("No 'Sign' label found.")
    else:       
        for shape in adjusted_shapes:
            if shape['label'] in ['Wai', 'Nei']:
                combined_array.extend(shape['points'])
                
    AN1 = create_convex_hull(np.array(combined_array))
    maskAN1 = create_mask_from_polygon(AN1, size) 
    d_array = caldc(maskB, maskAN1, k)
    rc, pards = calrc(d_array, epsilon=epsilon)
    data['shapes'] = adjusted_shapes
    return rc, pards, data, AN1

def main():
    parser = argparse.ArgumentParser(description="Registration coordinate update.")       
    parser.add_argument("-j", "--jsonFile", type=str, default='Adjusted_Output.json',
                        help="The JSON file with coordinates (default: Adjusted_Output.json)")
    parser.add_argument("-o", "--outputJsonFile", type=str, default='AllRegistrationSchemes.json',
                        help="The output JSON file with adjusted coordinates (default: AllRegistrationSchemes.json)")
    parser.add_argument("-cp", "--cellPointsFile", type=str, default='cellPoints.txt',
                        help="The file with cell points (default: cellPoints.txt)")                        
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
    parser.add_argument("-S", "--signal", action="store_true", help="Add Signal field if set") 
                            
    args = parser.parse_args()
    jsonFile = args.jsonFile
    outputJsonFile = args.outputJsonFile
    cellPointsFile = args.cellPointsFile    
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
    rotatetheta = np.max(thetas) 
    B = np.loadtxt(cellPointsFile, delimiter='\t', dtype=float)
    B = create_convex_hull(B)
    B = B.copy()

    maskB = create_mask_from_polygon(B, size)
    centerB = find_extrema_points(maskB)
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
    maskAN = create_mask_from_polygon(AN, size)
    centerA = find_extrema_points(maskAN) 

    adjusted_shapes = []
    for shape in data['shapes']:
        points = [rotate_points([x, y], centerA, transform=transform) for x, y in shape['points']]
        points = [rotate_point([x, y], rotatetheta, rotation_center=centerA) for x, y in points]
        points = translate_points(points, centerB, centerA)
        shape['points'] = points
        adjusted_shapes.append(shape)
    data['shapes'] = adjusted_shapes

    combined_array = []
    if args.signal:
        sign_points = None
        for shape in adjusted_shapes:
            if shape['label'] == 'Sign':
                sign_points = shape['points']
                break

        if sign_points is not None:
            combined_array = sign_points
        else:
            print("No 'Sign' label found.")
    else:       
        for shape in adjusted_shapes:
            if shape['label'] in ['Wai', 'Nei']:
                combined_array.extend(shape['points'])    

    AN = create_convex_hull(np.array(combined_array))
    AN = AN.copy()
    maskAN = create_mask_from_polygon(AN, size)
    centerA = find_extrema_points(maskAN) 
    
    convex_pointsA = get_convex_points(maskAN)
    dis_A1 = [calculate_distance(point, centerA) for point in convex_pointsA]
    mean_dis_A = np.mean(dis_A1)
    convex_pointsB = get_convex_points(maskB)
    dis_B1 = [calculate_distance(point, centerB) for point in convex_pointsB]
    mean_dis_B = np.mean(dis_B1) 
    dis0 = mean_dis_B / mean_dis_A
    data_to_store = []  
    dis_array = []  
    shapes = data['shapes']
    ratios1 = {
        'D': dis0, 
        'B': dis0, 
        'C': dis0,  
        'A': dis0  
    }
    
    dis_array.append(ratios1)
    results = Parallel(n_jobs=2)(delayed(process_theta_dis)(ratios, shapes, centerA, centerB, size, k, epsilon, maskB, args.signal) for ratios in dis_array)
    
    RC, pards, ratios, adjusted_shapesAll, Aall = zip(*results)
    filtered_RC = [r for r in RC if r is not None]
    filtered_pards = [p for r, p in zip(RC, pards) if r is not None]
    filtered_ratios = [r_ for r, r_ in zip(RC, ratios) if r is not None]
    filtered_adjusted_shapesAll = [a for r, a in zip(RC, adjusted_shapesAll) if r is not None]
    filtered_Aall = [a for r, a in zip(RC, Aall) if r is not None]

    for r, pas, ratio_dict, shapestmp, A in zip(filtered_RC, filtered_pards, filtered_ratios, filtered_adjusted_shapesAll, filtered_Aall):
        data['shapes'] = shapestmp
        tmp = {
            'RC': r,
            'PARDS': pas,
            'THETA': rotatetheta,
            'RATIOS': ratio_dict,
            'TRAN': transform,
            'JSON': deepcopy(data),
            'MASKPOINT': np.array(A, dtype=np.float32).tolist(),
            'CENT': np.array(centerB, dtype=np.float32).tolist()
        }            
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
            rc, pards, datatmp, Ad = process_dis(item_copy['JSON'], item_copy['CENT'], size, k, epsilon, ratiosUpdata, maskB, thetaIN, args.signal)
            tmpu = {
                'RC': rc,
                'PARDS': pards,
                'THETA': rotatetheta,
                'THETAIN': thetaIN,
                'RATIOS': ratio_dict,
                'RATIOSUpdate': ratiosUpdata,
                'TRAN': transform,
                'JSON': datatmp, 
                'MASKPOINT': np.array(Ad, dtype=np.float32).tolist(),
                'CENT': np.array(centerB, dtype=np.float32).tolist()
            }            
            data_to_storeNew.append(tmpu)

    data_to_storeF = data_to_store + data_to_storeNew
    json_data = json.dumps(data_to_storeF, indent=4) 

    with open(outputJsonFile, 'w') as json_file:
        json_file.write(json_data)

    print(f'Data has been written to ({outputJsonFile})')


if __name__ == "__main__":
    main()

