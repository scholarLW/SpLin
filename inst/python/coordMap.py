import argparse
import json
from adsco import adjust_shape_coordinates as asc

def main():
    parser = argparse.ArgumentParser(description="Coordinate Mapping Algorithm.")
    parser.add_argument("-j", "--jsonFile", type=str, default='Output.json',
                        help="The JSON file with coordinates (default: Output.json)")
    parser.add_argument("-o", "--outputJsonFile", type=str, default='Adjusted_Output.json',
                        help="The output JSON file with adjusted coordinates (default: Adjusted_Output.json)")
    parser.add_argument("-nw", "--newWidth", type=int, default=1170,
                        help="New width for the resized image (default: 1170)")
    parser.add_argument("-nh", "--newHeight", type=int, default=1108,
                        help="New height for the resized image (default: 1108)")
    parser.add_argument("-ow", "--oldWidth", type=int, default=23593,
                        help="Old width for the original image (default: 23593)")
    parser.add_argument("-oh", "--oldHeight", type=int, default=22342,
                        help="Old height for the original image (default: 22342)")
    args = parser.parse_args()

    with open(args.jsonFile, 'r') as file:
        data = json.load(file)

    adjusted_shapes = asc(data['shapes'], args.newWidth, args.newHeight, args.oldWidth, args.oldHeight)
    data['shapes'] = adjusted_shapes

    with open(args.outputJsonFile, 'w') as file:
        json.dump(data, file, indent=4)

    print(f"The coordinates have been updated and saved to{args.outputJsonFile}")

if __name__ == "__main__":
    main()