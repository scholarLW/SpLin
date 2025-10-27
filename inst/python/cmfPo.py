import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import cv2
import numpy as np
from matplotlib.patches import Polygon

def create_mask_from_polygon(polygon, size, transform=None):
    height, width = size[1], size[0]
    mask = np.zeros((height, width), dtype=np.uint8)
    cv2.fillPoly(mask, [polygon.astype(np.int32)], 255)

    if transform == 'horizontal':
        mask = cv2.flip(mask, 1)
    elif transform == 'vertical':
        mask = cv2.flip(mask, 0)

    return mask

def create_mask_from_image(image_path, transform=None):
    mask_image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

    if mask_image is None:
        raise ValueError(f"Unable to read the image file:{image_path}")

    _, mask = cv2.threshold(mask_image, 127, 255, cv2.THRESH_BINARY)

    if transform == 'horizontal':
        mask = cv2.flip(mask, 1)
    elif transform == 'vertical':
        mask = cv2.flip(mask, 0)

    return mask



 