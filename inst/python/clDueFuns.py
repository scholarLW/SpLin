import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import cv2
import math
import numpy as np


def estimate_transformation(src_pts, dst_pts):
    tfm = cv2.getPerspectiveTransform(np.array(src_pts, dtype=np.float32),
                                      np.array(dst_pts, dtype=np.float32))
    return tfm

def apply_transformation(points, transformation):
    points = np.array([points], dtype=np.float32)
    transformed_points = cv2.perspectiveTransform(points, transformation)
    return transformed_points[0]

def calculate_centroid(vertices):
    num_vertices = len(vertices)
    sum_x = sum(vertex[0] for vertex in vertices)
    sum_y = sum(vertex[1] for vertex in vertices)

    centroid_x = sum_x / num_vertices
    centroid_y = sum_y / num_vertices

    return [centroid_x, centroid_y]

def calculate_min(vertices):
    min_x = min(vertex[0] for vertex in vertices)
    return min_x

def calculate_intersection_point(A, B, C):
    x_A, y_A = A
    x_B, y_B = B
    x_C, y_C = C

    if x_A == x_B:
        x_O = x_A
        y_O = y_C
    else:
        if y_A == y_B:
            x_O = x_C
            y_O = y_A
        else:
            m_AB = (y_B - y_A) / (x_B - x_A)
            m_per = -1 / m_AB
            x_O = (y_A * x_B - y_B * x_A - (x_B - x_A) * (y_C - m_per * x_C)) / ((x_B - x_A) * m_per + y_A - y_B)
            y_O = m_per * x_O + y_C - m_per * x_C

    return (x_O, y_O)

def point_side_of_line(A, C, O):
    vector_AC = (C[0] - A[0], C[1] - A[1])
    vector_AO = (O[0] - A[0], O[1] - A[1])

    cross_product = vector_AC[0] * vector_AO[1] - vector_AC[1] * vector_AO[0]

    if cross_product > 0:
        return "left"
    elif cross_product < 0:
        return "right"
    else:
        return "on_line"

def calculate_angle_BINI(x_O, y_O, x_P, y_P, type=0):
    length_CO = math.sqrt((x_P - x_O) ** 2 + (y_P - y_O) ** 2)
    cos_value = (y_P - y_O) / length_CO
    angle_radians = math.acos(cos_value)
    if type == 0:
        final_angle = math.pi + angle_radians
    else:
        if type == 1:
            final_angle = angle_radians
        else:
            if type == 2:
                final_angle = 2 * math.pi - angle_radians
            else:
                final_angle = math.pi - angle_radians

    return final_angle

def calculate_thetha_BINI(C, A, P, O):
    pointsside = point_side_of_line(A, C, O)
    if pointsside == 'on_line':
        print("Error!")
        sys.exit(1)

    if (C[1] <= A[1] and C[0] <= A[0] and pointsside == 'right') or (C[1] <= A[1] and C[0] > A[0] and pointsside == 'right'):
        thetha = calculate_angle_BINI(O[0], O[1], P[0], P[1], type=0)
    else:
        if (C[1] > A[1] and C[0] > A[0] and pointsside == 'left') or (C[1] > A[1] and C[0] <= A[0] and pointsside == 'left'):
            thetha = calculate_angle_BINI(O[0], O[1], P[0], P[1], type=1)
        else:
            if (C[1] <= A[1] and C[0] <= A[0] and pointsside == 'left') or (C[1] <= A[1] and C[0] > A[0] and pointsside == 'left'):
                thetha = calculate_angle_BINI(O[0], O[1], P[0], P[1], type=2)
            else:
                thetha = calculate_angle_BINI(O[0], O[1], P[0], P[1], type=3)

    return thetha

def translate_points_per(points, A, AF):
    A = np.array([float(item) for item in A], dtype=float)
    AF = np.array([float(item) for item in AF], dtype=float)
    points = np.array([float(item) for item in points], dtype=float)

    dx = A[0] - AF[0]
    dy = A[1] - AF[1]

    translated_points = [points[0] + dx, points[1] + dy]
    return translated_points

def point_to_line_distance(A, B, P):
    A = np.array([float(item) for item in A])
    B = np.array([float(item) for item in B])
    P = np.array([float(item) for item in P])
    distance = np.float32(abs((B[0] - A[0]) * P[1] + (A[1] - B[1]) * P[0] + B[1] * A[0] - A[1] * B[0]) / math.sqrt(
        (B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2))
    return distance

def calculate_angle(x_O, y_O, x_P, y_P, type=0):
    length_CO = math.sqrt((x_P - x_O) ** 2 + (y_P - y_O) ** 2)

    cos_value = (y_P - y_O) / length_CO

    angle_radians = math.acos(cos_value)

    if type == 0:
        final_angle = math.pi - angle_radians
    else:
        final_angle = math.pi + angle_radians

    return final_angle

def calculate_thetha(A, B, P, O):
    if (A[0] <= B[0] and A[1] >= B[1] and O[1] < P[1]) or (A[0] <= B[0] and A[1] <= B[1] and O[1] >= P[1]) or (A[0] >= B[0] and A[1] <= B[1] and O[1] < P[1]) or (A[0] >= B[0] and A[1] >= B[1] and O[1] < P[1]):
        thetha = calculate_angle(O[0], O[1], P[0], P[1], type=0)
    else:
        thetha = calculate_angle(O[0], O[1], P[0], P[1], type=1)

    return thetha

