import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import cv2
import numpy as np

def find_extrema_points(mask):
    _, thresh = cv2.threshold(mask, 254, 255, cv2.THRESH_BINARY)
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None, None, None

    M = cv2.moments(contours[0])
    if M["m00"] == 0:
        return None, None, None
    cX = np.float32(M["m10"] / M["m00"])
    cY = np.float32(M["m01"] / M["m00"])

    return [cX, cY]