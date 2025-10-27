import argparse
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = "20000000000"
import json
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from ssDRe import ssDNAReSize as SDR
from heReS import HEReSize as HER
from sigNal import SigNSize as SGN
from PIL import Image

def save_info_to_file(filename, original_width, original_height, new_width, new_height):
    image_info = {
        "scaled_dimensions": {
            "width": new_width,
            "height": new_height
        },
        "original_dimensions": {
            "width": original_width,
            "height": original_height
        }
    }

    with open(filename, 'w') as file:
        json.dump(image_info, file, indent=4)
       

def main():
    parser = argparse.ArgumentParser(description="Resize an image and save its dimensions.")
    parser.add_argument("-I", "--inputFile", type=str, default='./CTR95_HE.tiff',
                        help="The original image file (default: ./CTR95_HE.tiff)")
    parser.add_argument("-O", "--outputFile", type=str, default=None,
                        help="Save the resized image file (default: None)")
    parser.add_argument("-d", "--dpi", type=int, default=300,
                        help="DPI for the resized image (default: 300)")
    parser.add_argument("-t", "--imageType", type=str, choices=['ssDNA', 'HE', 'Signal'], default='HE',
                        help="Type of the original image (default: Signal)")
    args = parser.parse_args()

    Image.MAX_IMAGE_PIXELS = None 
    
    if args.outputFile is not None:
        if not os.path.dirname(args.outputFile):
            args.outputFile = os.path.join(os.getcwd(), args.outputFile)

        directory = os.path.dirname(args.outputFile)
    else:
        directory = os.path.dirname(args.inputFile)
        
    if not directory:
        directory = os.getcwd()
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

    if args.imageType == 'HE':
        original_width, original_height, new_width, new_height = HER(args.inputFile, args.outputFile, args.dpi)
    if args.imageType == 'ssDNA':
        original_width, original_height, new_width, new_height = SDR(args.inputFile, args.outputFile, args.dpi)
    if args.imageType == 'Signal':
        original_width, original_height, new_width, new_height = SGN(args.inputFile)

    print(f"The width and height of the rescaled image are ({new_width}, {new_height}), corresponding to the original image's width and height of ({original_width}, {original_height}).")

    outputTxtFile = os.path.join(directory, 'image_dimensions.json')
    save_info_to_file(outputTxtFile, original_width, original_height, new_width, new_height)
    
    print(f"The information has been saved to {outputTxtFile}.")


if __name__ == "__main__":
    main()