import argparse
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import sys
import csv
import cv2
import math
import numpy as np
from tranDueFuns import rotate_point, calculate_distance
from clDueFuns import estimate_transformation, apply_transformation, calculate_centroid, calculate_min, calculate_intersection_point
from clDueFuns import point_side_of_line, calculate_angle_BINI, calculate_thetha_BINI, translate_points_per, point_to_line_distance, calculate_angle, calculate_thetha


def write_dict_list_to_txt(data, filename):
    with open(filename, 'w', encoding='utf-8') as file:
        headers = "Cell\tx\ty\tbin\txnew\tynew\n"
        file.write(headers)

        for row in data:
            row_data = "\t".join(str(row[key]) for key in ["Cell", "x", "y", "bin", "xnew", "ynew"])
            file.write(row_data + "\n")


def main():
    parser = argparse.ArgumentParser(description="Linearization and straightening of spatial coordinates.")
    parser.add_argument("-mdp", "--metaDtaPointsFile", type=str, default='metaDtaPoints.txt',
                        help="The metadata points file (default: metaDtaPoints.txt)")
    parser.add_argument("-o", "--outputFile", type=str, default='metaDtaPointsUpdate.txt',
                        help="The output file for updated metadata points (default: metaDtaPointsUpdate.txt)")
    parser.add_argument("-bfp", "--BinFPointsFile", type=str, default=None,
                        help="The binF points file (default: None)")
    parser.add_argument("-bip", "--BinIPointsFile", type=str, default='BinIPoints.txt',
                        help="The binI points file (default: BinIPoints.txt)")    

    args = parser.parse_args()
    metaDtaPointsFile = args.metaDtaPointsFile
    outputFile = args.outputFile
    BinFPointsFile = args.BinFPointsFile
    BinIPointsFile = args.BinIPointsFile
    
    scale = False
    
    BinIPoints = []

    with open(BinIPointsFile, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t') 
        for row in reader:
            BinIPoints.append(row)

    metaDtaPoints = []
    with open(metaDtaPointsFile, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            metaDtaPoints.append(row)

    HighNUM = []
    BinIPointsNEW = []
    adjusted_metaDtaPoints = []
    irow = 1
    BNEWtmp = []
    for row in BinIPoints:
        Dold = BNEWtmp

        src_pts = [np.float32(row[key]) for key in ['Cx', 'Cy', 'Ax', 'Ay', 'Bx', 'By', 'Dx', 'Dy']]
        dst_pts = [np.float32(row[key]) for key in ['Cx1', 'Cy1', 'Ax1', 'Ay1']]
        centroidO = np.array([np.float32(row[key]) for key in ['Ox', 'Oy']])

        src_points = np.array([src_pts[i:i + 2] for i in range(0, len(src_pts), 2)], dtype=np.float32)
        dst_points = np.array([dst_pts[i:i + 2] for i in range(0, len(dst_pts), 2)], dtype=np.float32)
  
        P = calculate_intersection_point(src_points[0], src_points[1], centroidO)

        theta = calculate_thetha_BINI(src_points[0], src_points[1], P, centroidO)

        PNEW = rotate_point(P, theta, rotation_center=centroidO)
        CNEW = rotate_point(src_points[0], theta, rotation_center=centroidO)
        BNEW = rotate_point(src_points[2], theta, rotation_center=centroidO)
        DNEW = rotate_point(src_points[3], theta, rotation_center=centroidO)
        BNEWtmp = [BNEW[0], BNEW[1]]
        DNEWtmp = [DNEW[0], DNEW[1]]

        pointsside = point_side_of_line(src_points[1], src_points[0], centroidO)
        if pointsside == 'on_line':
            print("Error!")
            sys.exit(1)

        if (pointsside == 'left'):
            BNEWtmp[1] = 2 * centroidO[1] - BNEW[1]
            DNEWtmp[1] = 2 * centroidO[1] - DNEW[1]

        if PNEW[0] < CNEW[0]:
            centroidO_New_X = [dst_points[0][0] - calculate_distance(CNEW, PNEW), calculate_distance(centroidO, P)]
        else:
            centroidO_New_X = [dst_points[0][0] + calculate_distance(CNEW, PNEW), calculate_distance(centroidO, P)]

        BNEWtmp = translate_points_per(np.float32(BNEWtmp), centroidO_New_X, centroidO)
        DNEWtmp = translate_points_per(np.float32(DNEWtmp), centroidO_New_X, centroidO)

        ODIS = point_to_line_distance(src_points[0], src_points[1], centroidO)
        HighNUM.append(ODIS)

        for rowm in metaDtaPoints:
            if rowm['bin'] == row['ID']:
                x = rowm['x']
                y = rowm['y']
                points = [x, y]
                adjusted_points = rotate_point(points, theta, rotation_center=centroidO)
                new_x, new_y = np.float32(adjusted_points)

                if (pointsside == 'left'):
                    new_y = 2 * centroidO[1] - new_y

                adjusted_points = translate_points_per([new_x, new_y], centroidO_New_X, centroidO)
                new_x, new_y = np.float32(adjusted_points)
                rowm['xnew'] = new_x
                rowm['ynew'] = new_y

                adjusted_metaDtaPoints.append(rowm)

        if irow == 1:
            Dold = DNEWtmp

        BinIPointsNEW.append(
            {'ID': row['ID'], 'Cx': dst_pts[0], 'Cy': dst_pts[1], 'Ax': dst_pts[2], 'Ay': dst_pts[3], 'Bx': BNEWtmp[0],
             'By': BNEWtmp[1], 'Dx': DNEWtmp[0], 'Dy': DNEWtmp[1], 'Cx1': dst_pts[0], 'Cy1': dst_pts[1],
             'Ax1': dst_pts[2], 'Ay1': dst_pts[3], 'Bx1': BNEWtmp[0], 'By1': BNEWtmp[1], 'Dx1': Dold[0],
             'Dy1': Dold[1]})

        irow = irow + 1

    if scale:
        adjusted_metaDtaPoints_up = []
        for row in BinIPointsNEW:
            src_pts = [np.float32(row[key]) for key in ['Cx', 'Cy', 'Ax', 'Ay', 'Bx', 'By', 'Dx', 'Dy']]
            dst_pts = [np.float32(row[key]) for key in ['Cx1', 'Cy1', 'Ax1', 'Ay1', 'Bx1', 'By1', 'Dx1', 'Dy1']]

            src_points = np.array([src_pts[i:i + 2] for i in range(0, len(src_pts), 2)], dtype=np.float32)
            dst_points = np.array([dst_pts[i:i + 2] for i in range(0, len(dst_pts), 2)], dtype=np.float32)

            transformation = estimate_transformation(src_points, dst_points)

            for rowm in adjusted_metaDtaPoints:
                if rowm['bin'] == row['ID']:
                    x = rowm['xnew']
                    y = rowm['ynew']
                    points = [x, y]
                    adjusted_points = apply_transformation([points], transformation)
                    new_x, new_y = adjusted_points[0]
                    rowm['xnew'] = new_x
                    rowm['ynew'] = new_y
                    adjusted_metaDtaPoints_up.append(rowm)
        adjusted_metaDtaPoints = adjusted_metaDtaPoints_up

    if BinFPointsFile is not None and os.path.exists(BinFPointsFile):
        BinFPoints = np.genfromtxt(BinFPointsFile, delimiter='\t', missing_values='NA', filling_values=np.nan, skip_header=1)

        if BinFPoints.size > 0:
            nbinI = len(BinIPointsNEW)
            binIA = [[row['Ax1'], row['Ay1']] for row in [BinIPointsNEW[nbinI - 1]]][0]
            binIB = [[row['Bx1'], row['By1']] for row in [BinIPointsNEW[nbinI - 1]]][0]
            centroidO = calculate_centroid(BinFPoints)

            P = [calculate_min(BinFPoints), centroidO[1]]

            centroidO_New_X = [np.float32(max(binIA[0], binIB[0])) + 2 * calculate_distance(centroidO, P), 0]

            for rowm in metaDtaPoints:
                if rowm['bin'] == 'binF':
                    x = rowm['x']
                    y = rowm['y']
                    points = [x, y]

                    adjusted_points = translate_points_per(points, centroidO_New_X, centroidO)
                    new_x, new_y = np.float32(adjusted_points)
                    rowm['xnew'] = new_x
                    rowm['ynew'] = new_y

                    adjusted_metaDtaPoints.append(rowm)

    write_dict_list_to_txt(adjusted_metaDtaPoints, outputFile)
    print("The file has been written.")


if __name__ == "__main__":
    main()

