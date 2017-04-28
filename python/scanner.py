import numpy as np
import scipy
import matplotlib.pylab as plt
import cornerDetector

def scanner(image, displayFigures, axs):
  # Local Variables: y2, target_points, lengthY, lengthX, target, displayFigures, H, image, v, y, x2, points, u, x, y1, x3, y3, x1, y4, x4
  # Function calls: scanner, homography_solve, uint8, homography_transform, cornerDetector, zeros, round, size
  # Find corners in original image, and transform using homography
  # The input arguments are:
  #    image = original image to be transformed
  #    angle = angle between camera and normal surface of the plane
  # The output argument is:
  #   target = transformed image

  #Find corners
  points = cornerDetector.cornerDetector(image, displayFigures, axs)
  target = np.zeros((3, 3, 3), np.uint8)
  return target
'''
  # Constrain the corner
  lengthX = round(points[1,2] - points[1,3]) + 1
  lengthY = round(1.3*lengthX) #round(1.1*(points(2,2) - points(2,1))) + 1;
  x4 = 1
  y4 = 1
  x1 = x4 + lengthX
  y1 = y4
  x3 = x4
  y3 = y4 + lengthY
  x2 = x1
  y2 = y3
  target_points = np.array([[x1, x2, x3, x4],
                            [y1, y2, y3, y4]])
                   
  H = homography_solve(target_points, points);
'''
#add target back here
'''
  for x in range(0,lengthX)
    for y in range(0, lengthY)
      u = np.array([[x],[y]])
      v = homography_transform(u, H);
      target[y,x] = image[round(v[2]),round(v[1])];
'''
