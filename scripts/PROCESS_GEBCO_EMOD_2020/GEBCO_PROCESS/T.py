import iris
import numpy as np
import argparse
from pathlib import Path
import sys


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
    sys.exit("Sorry all the arguments are required") 


print (args.AMM15_PATH[0])
print("\n----------------------------------------------------\n")
print( "Thanks, you have chosen: \n ")
print( "      AMM15 grid file as {}\n".format (args.AMM15_PATH[0]))
if Path(args.AMM15_PATH[0]).is_file():
       print(" and the file {} exists.".format (args.AMM15_PATH[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.AMM15_PATH[0])) 

print("\n     the INPUT directory for bathy as {}\n".format(  args.IN_DIR[0] ))
if (Path(args.IN_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.IN_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.IN_DIR[0])) 
print("\n     the OUTPUT directory for processed bathy as {}\n".format(  args.OUT_DIR[0] ))
if (Path(args.OUT_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.OUT_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.OUT_DIR[0])) 

#       and output for process bathy as {}.".format(args.AMM15_PATH[0], args.IN_DIR[0], args.OUT_DIR[0]) ) 
print("\n----------------------------------------------------\n")

AMM15_cube = iris.load(args.AMM15_PATH)[0]

print(AMM15_cube.coord_dims)
print(AMM15_cube.coords)
print(AMM15_cube.coord_system)
X_LAT= np.array((AMM15_cube.coord('grid_latitude').points[:]))
X_LON= np.array((AMM15_cube.coord('grid_longitude').points[:]))








