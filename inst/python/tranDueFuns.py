import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import cv2
import math
import numpy as np

def translate_points(points, A, AF):
    dx = A[0] - AF[0]
    dy = A[1] - AF[1]
    translated_points = [[x + dx, y + dy] for x, y in points]
    return translated_points
    
def due_points(point, center_from, center_to):
    x, y = point
    x_from, y_from = center_from
    x_to, y_to = center_to

    offset_x = x - x_from
    offset_y = y - y_from

    new_x = x_to + offset_x
    new_y = y_to + offset_y

    return [new_x, new_y]

def find_convex_hull(contour):
    hull = cv2.convexHull(contour, returnPoints=False)
    points = cv2.convexHull(contour, hull)
    return points

def get_convex_points(mask):
    _, binary_mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
    contours, _ = cv2.findContours(binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    points_list = []
    for contour in contours:
        for point in contour:
            x, y = point[0]
            points_list.append([x, y])

    return points_list
    
def calculate_distance(point1, point2):
    point1 = np.array([float(item) for item in point1], dtype=float)
    point2 = np.array([float(item) for item in point2], dtype=float)
    x1, y1 = point1
    x2, y2 = point2
    distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    return distance
   
def calculate_rotation_angle(op, op_prime):
    x1, y1 = op
    x2, y2 = op_prime
    dot_product = x1 * x2 + y1 * y2
    mod_op = math.sqrt(x1 ** 2 + y1 ** 2)
    mod_op_prime = math.sqrt(x2 ** 2 + y2 ** 2)
    cos_theta = dot_product / (mod_op * mod_op_prime)
    cos_theta = max(min(cos_theta, 1), -1)
    theta_radians = math.acos(cos_theta)
    return theta_radians   
    
def rotate_point(point, theta_rad, rotation_center=(0, 0)):
    rotation_matrix = np.array([
        [np.cos(theta_rad), -np.sin(theta_rad)],
        [np.sin(theta_rad), np.cos(theta_rad)]
    ])
    point = np.array(point, dtype=float)
    rotation_center = np.array(rotation_center, dtype=float)
    point_relative = point - rotation_center
    rotated_point_relative = np.dot(rotation_matrix, point_relative)
    rotated_point = rotated_point_relative + rotation_center

    return tuple(rotated_point)

def rotate_points(points, centroidO, transform=None):
    x = points[0]
    y = points[1]
    if transform == 'horizontal':
        x = 2 * centroidO[0] - x
    elif transform == 'vertical':
        y = 2 * centroidO[1] - y

    return [x, y]
    
def scale_polygon_points(points, scale_factor, centroid):
    center = np.array(centroid, dtype=float)
    points = np.array(points, dtype=float)
    scaled_points = (points - center) * scale_factor + center

    return scaled_points.tolist()
 
def scale_polygon_points_NEW(points, ratios, centroid):
    center = np.array(centroid, dtype=float)
    points = np.array(points, dtype=float)
    A = ratios['A']
    B = ratios['B']
    C = ratios['C']
    D = ratios['D']
    scaled_points = []
    for point in points:
        x, y = point
        scale_factor_x = B if x > center[0] else D
        scale_factor_y = A if y > center[1] else C
        scaled_point = (point - center) * [scale_factor_x, scale_factor_y] + center
        scaled_points.append(scaled_point.tolist())
    scaled_points = np.array(scaled_points, dtype=float).tolist()
    
    return scaled_points    