import argparse
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import json
import cv2
import numpy as np

def convert_coordinate(original_width, original_height, new_width, new_height, x_resized, y_resized):
    x_original = round(x_resized * original_width / new_width, 13)
    y_original = original_height - round(y_resized * original_height / new_height, 13)
    return [x_original, y_original]

def adjust_shape_coordinates(shapes, neww, newh, oldw, oldh):
    adjusted_shapes = []
    for row in shapes:
        x = float(row['X'])   
        y = float(row['Y'])  
        adjusted_points = convert_coordinate(oldw, oldh, neww, newh, x, y)
        row['X'] = adjusted_points[0]
        row['Y'] = adjusted_points[1]
        adjusted_shapes.append(row)
    return adjusted_shapes


def main():
    parser = argparse.ArgumentParser(description="Extract pixels and colors from images.")
    parser.add_argument("-j", "--jsonFile", type=str, default='DSS12_ssDNA.json',
                        help="The JSON file with coordinates (default: DSS12_ssDNA.json)")
    parser.add_argument("-i", "--imageFile", type=str, default='DSS12_ssDNA.png',
                        help="The scaled image file (default: DSS12_ssDNA.png)")            
    parser.add_argument("-o", "--outputCoordinatesFile", type=str, default='Adjusted_output_coordinates_and_colors.txt',
                        help="The adjusted output file for coordinates and colors (default: Adjusted_output_coordinates_and_colors.txt)")
    parser.add_argument("-ww", "--newWidth", type=int, default=1170,
                        help="The new width of the scaled image (default: 1170)")
    parser.add_argument("-wh", "--newHeight", type=int, default=1108,
                        help="The new height of the scaled image (default: 1108)")
    parser.add_argument("-ow", "--oldWidth", type=int, default=23593,
                        help="The old width of the original image (default: 23593)")
    parser.add_argument("-oh", "--oldHeight", type=int, default=22342,
                        help="The old height of the original image (default: 22342)")
    
    args = parser.parse_args()
    jsonFile = args.jsonFile
    imageFile = args.imageFile
    outputCoordinatesFile = args.outputCoordinatesFile
    neww = args.newWidth
    newh = args.newHeight
    oldw = args.oldWidth
    oldh = args.oldHeight

    directory = os.path.dirname(outputCoordinatesFile)

    if not directory:
        directory = os.getcwd()
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

    try:
        with open(jsonFile, 'r') as file:
            data = json.load(file)
    except FileNotFoundError:
        print(f"Error: File {jsonFile} not found.")
        return
    except json.JSONDecodeError:
        print(f"Error: Cannot decode {jsonFile} as JSON.")
        return

    if len(data['shapes']) > 1:
        points_wai = None
        points_nei = None

        for shape in data['shapes']:
            if shape['label'] == 'Wai':
                points_wai = shape['points']
            elif shape['label'] == 'Nei':
                points_nei = shape['points'][::-1]  

        if points_wai is None or points_nei is None:
            print("Error: Missing 'Wai' or 'Nei' label in the JSON file.")
            return

        if points_wai[-1] == points_nei[0]:
            points = points_wai + points_nei[1:]
        else:
            points = points_wai + points_nei

        points = np.array(points, dtype=np.int32)
    else:
        points = np.array(data['shapes'][0]['points'], dtype=np.int32)

    if not np.array_equal(points[0], points[-1]):
        points = np.vstack([points, points[0]])

    image = cv2.imread(imageFile)
    if image is None:
        print(f"Error: File {imageFile} not found or could not be read.")
        return

    mask = np.zeros(image.shape[:2], dtype=np.uint8)
    cv2.fillPoly(mask, [points], 255)

    masked_image = cv2.bitwise_and(image, image, mask=mask)

    cv2.imwrite(os.path.join(directory, 'masked_image.png'), masked_image)

    pixels_list = []
    for y in range(masked_image.shape[0]):
        for x in range(masked_image.shape[1]):
            if masked_image[y, x].any(): 
                color = tuple(masked_image[y, x])
                pixel_info = {
                    "X": x,
                    "Y": y,
                    "Color": f"#{color[2]:02x}{color[1]:02x}{color[0]:02x}"
                }
                pixels_list.append(pixel_info)

    adjusted_shapes = adjust_shape_coordinates(pixels_list, neww, newh, oldw, oldh)

    with open(outputCoordinatesFile, "w") as f:
        f.write("X\tY\tColor\n")
        for item in adjusted_shapes:
            line = f"{item['X']:.6f}\t{item['Y']:.6f}\t{item['Color']}\n"
            f.write(line)

    print(f"Coordinates and colors have been written to {outputCoordinatesFile}")

if __name__ == "__main__":
    main()
