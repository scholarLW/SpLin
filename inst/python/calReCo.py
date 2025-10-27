import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import numpy as np
import cv2
import matplotlib.pyplot as plt

def intersection(maskImage):
    image = cv2.imread(maskImage, cv2.IMREAD_GRAYSCALE)
    if image is None:
        print("Image not found or unable to read.")
        return None, None
    _, mask = cv2.threshold(image, 127, 255, cv2.THRESH_BINARY)
    moments = cv2.moments(mask)
    if moments["m00"] != 0:
        cx = int(moments["m10"] / moments["m00"])
        cy = int(moments["m01"] / moments["m00"])
    else:
        print("Mask is empty.")
        return None, mask


    x_coords = []
    y_coords = []
    for x in range(mask.shape[1]):
        if mask[cy, x] == 255:
            x_coords.append(x)
    for y in range(mask.shape[0]):
        if mask[y, cx] == 255:
            y_coords.append(y)

    x_min = min(x_coords) if x_coords else None
    x_max = max(x_coords) if x_coords else None
    y_min = min(y_coords) if y_coords else None
    y_max = max(y_coords) if y_coords else None

    extreme_points = {
        'D': [x_min, cy],
        'B': [x_max, cy],
        'C': [y_min, cx],
        'A': [y_max, cx]
    }

    return extreme_points

def calculate_deviation_coefficients(image1, image2, k):
    assert image1.shape == image2.shape, "The pixel sizes of the two images are inconsistent."
    size = image1.shape[0]*image1.shape[1]
    contours1, _ = cv2.findContours(image1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours2, _ = cv2.findContours(image2, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    deviation_coefficients = []

    for y in range(0, image1.shape[0], k):
        intersection_points1 = []
        intersection_points2 = []

        for contour1, contour2 in zip(contours1, contours2):
            for point1 in contour1:
                if point1[0][1] == y:
                    intersection_points1.append(point1[0])
            for point2 in contour2:
                if point2[0][1] == y:
                    intersection_points2.append(point2[0])

        if len(intersection_points1) > 0 and len(intersection_points2) > 0:
            pointA1 = intersection_points1[0]
            pointA2 = intersection_points2[0]
            if len(intersection_points1) > 1:
                pointB1 = max(intersection_points1, key=lambda x: x[0])
            else:
                pointB1 = pointA1
            if len(intersection_points2) > 1:
                pointB2 = max(intersection_points2, key=lambda x: x[0])
            else:
                pointB2 = pointA2

            dA = abs(pointA1[0]-pointA2[0])*y/size
            dB = abs(pointB1[0]-pointB2[0])*y/size
            d = (dA + dB) / 2
            deviation_coefficients.append(d)

    return deviation_coefficients

def calculate_registration_coefficient(distance_diffs, epsilon=0.1):
    try:
        d = np.array(distance_diffs)

        if d.size == 0:
            print("The input list of distance differences is empty.")
            return None

        d_transformed = np.log1p(d) 
        pards = sum(d_transformed)
        
        Q1 = np.percentile(d_transformed, 25)
        Q3 = np.percentile(d_transformed, 75)

        IQR = Q3 - Q1

        mu_prime = np.median(d_transformed)

        if IQR == 0:
            return 1.0, pards
        else:
            registration_coefficient = max(0, 1 - np.log(1 + (1-epsilon)*mu_prime/IQR))
            return registration_coefficient, pards

    except Exception as e:
        print(f"Error occurred while calculating the registration coefficient:{e}")
        return None
