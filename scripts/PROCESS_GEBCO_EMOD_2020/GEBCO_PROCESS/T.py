import iris
import numpy as np
import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description='Process inputs and output file paths')
parser.add_argument('-a','--AMM15_PATH', metavar='AMM15_cube_file', nargs=1,
                    help='File location of AMM15 rotated grid', required=True)
parser.add_argument('-i','--IN_DIR',metavar='IN_DIR',  nargs=1,
                    help='Path to source bathymertry files', required=True)
parser.add_argument( '-o','--OUT_DIR',metavar='OUT_DIR', nargs=1,
                    help='Path to output to ', required=True) 

args = parser.parse_args()
if not all([args.AMM15_PATH, args.IN_DIR, args.OUT_DIR]):
    print(" Sorry All Arguments are required")
    exit(1)


print (args.AMM15_PATH[0])
print("\n----------------------------------------------------\n")
print( ("Thanks, you have chosen: \n \
       AMM15 grid file as {}\n \
       INPUT directory for bathy as {}\n \
       and output for process bathy as {}.".format(args.AMM15_PATH[0], args.IN_DIR[0], args.OUT_DIR[0]) ) )
print("\n----------------------------------------------------\n")

AMM15_cube = iris.load(args.AMM15_PATH)[0]

print(AMM15_cube.coord_dims)
print(AMM15_cube.coords)
print(AMM15_cube.coord_system)
X_LAT= np.array((AMM15_cube.coord('grid_latitude').points[:]))
X_LON= np.array((AMM15_cube.coord('grid_longitude').points[:]))








