import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from checker import Checker
import numpy as np
def plot(polygons):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for p in polygons:
        points = p.points
        for i in range(3):
            for j in range(i + 1, 3):
                ax.plot(*[[points[i][k], points[j][k]] for k in range(3)], color='gray')
    return fig, ax

def paint_and_save(fname):
    checker = Checker(fname)
    polygons = checker.model
    fig, ax = plot(polygons)
    fig.savefig(f'{fname[:fname.index(".stl")]} model')    
    fig.show()
    
if __name__ == '__main__':
    if len(sys.argv) > 2:
        fname = ' '.join(sys.argv[1:])
    else:
        fname = sys.argv[1]
    paint_and_save(fname)