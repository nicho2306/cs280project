import numpy as np
from scipy import ndimage as ndi
import matplotlib.pylab as plt
import matplotlib as mpl
from skimage import feature
from skimage import transform
import math
    
def rgb2gray(rgb):
  return np.dot(rgb[...,:3], [0.299, 0.587, 0.114])

def dist((x1,y1), (x2,y2)):
  return math.sqrt((x1-x2)**2 + (y1-y2)**2)

def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

# Return true if line segments AB and CD intersect
def intersect(A,B,C,D):
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def perp( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def seg_intersect(a1,a2, b1,b2) :
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b1

def plot_grayscale(axes, data):
    axes.imshow(data, norm=mpl.colors.Normalize(np.min(data), np.max(data)), cmap=plt.get_cmap('gray'))

def cornerDetector(image, displayFigures, axs):
  # Find corners in original image
  # The input arguments are:
  #    image = original image to be filtered
  #    displayCorners = whether to display image with plotted 4 corners
  #    displayHough   = whether to display image with lines from Hough transform 
  # The output argument is:
  #    points = array of (x,y) coordinates of the 4 corners detected
  
  ## convert to grayscale, lowpass filter, edge filter
  I = rgb2gray(image);
  I = ndi.filters.gaussian_filter(I, 10)
  BW = feature.canny(I, sigma=0.9)
  
  ## Hough transform
  h, theta, d = transform.hough_line(BW)
  lines = transform.probabilistic_hough_line(BW, line_length=110, line_gap=5)
  #_, angles, dists = hough_line_peaks(H,T,R, num_peaks=10, threshold=ceil(0.3*max(H(:))), min_disance=10000000000)
  #P = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
  #LINES = houghlines(BW,T,R,P, 'FillGap',10000000000,'MinLength',1);

  #merge similar angle and location lines
  slopeThresh=0.1
  distThresh = max(image.shape[0], image.shape[1])/3
  maxslope=0
  newLines = []
  for line in lines:
    newLines.append(line)
    dy=1.0*((line[1][1] - line[0][1]))
    dx=1.0*((line[1][0] - line[0][0]))
    if dx == 0.:
      dx = 1.0
    if(np.abs(dy/dx) > maxslope):
      maxslope = np.abs(dy/dx)

  for lineLoop in range(3):
    print len(newLines)
    '''
    plt.figure()
    for line in newLines[0:10]:
      ys = np.array([line[0][0],line[1][0]])
      xs = np.array([line[0][1],line[1][1]])
      plt.plot(ys, xs)
    plt.show()
    '''
    for line0 in lines:
      dy0=1.0*((line0[1][1] - line0[0][1]))
      dx0=1.0*((line0[1][0] - line0[0][0]))
      if dx0 == 0.:
        dx0 = 1.0
      #print dy0/dx0
      for line1 in newLines:
        if(line1 is line0):
          continue
        if(line0 not in newLines):
          continue
        if(line1 not in newLines):
          continue
        #print "intersect:" + str(intersect(np.asarray(line0[0]),np.asarray(line0[1]),np.asarray(line1[0]),np.asarray(line1[1])))
        dy1=1.0*((line1[1][1] - line1[0][1]))
        dx1=1.0*((line1[1][0] - line1[0][0]))
        if dx1 == 0.:
          dx1 = 1.0
        print "slopes:"+str(dy1/dx1/maxslope)+":0:"+str(dy0/dx0/maxslope)
        if(np.abs(dy1/dx1 - dy0/dx0)/maxslope < slopeThresh):
          d00 = dist(line1[0],line0[0])
          d01 = dist(line1[0],line0[1])
          d10 = dist(line1[1],line0[0])
          d11 = dist(line1[1],line0[1])
          mind = min(d00,d01,d10,d11)
          print "dist:"+str(mind)
          print line0
          print line1
          #check which endpoint is close\
          if(mind < distThresh):
            newLines.remove(line1)
            newLines.remove(line0)
            if(mind is d00):
              print "app 0,0"
              newLines.append((line1[1], line0[1]))
            elif(mind is d01):
              print "app 0,1"
              newLines.append((line1[1], line0[0]))
            elif(mind is d10):
              print "app 1,0"
              newLines.append((line1[0], line0[1]))
            elif(mind is d11):
              print "app 1,1"
              newLines.append((line1[0], line0[0]))
          else:
            pass #parallel lines that aren't close
        else: 
          pass # not the same slop so we don't care if they are close
  
  print lines
  print newLines
  ## Plot Hough lines & corners
  
  if displayFigures:
      plot_grayscale(axs[1], BW)
      for line in lines:
        ys = np.array([line[0][0],line[1][0]])
        xs = np.array([line[0][1],line[1][1]])
        axs[1].plot(ys, xs)
        axs[2].plot(ys, xs)
      for line in newLines:
        ys = np.array([line[0][0],line[1][0]])
        xs = np.array([line[0][1],line[1][1]])
        axs[3].plot(ys, xs)
  
  '''
  ## Find equation of lines
  max_area = 0;
  max_points = zeros(2,4);
  
  # random sample ~210 combinations of 10 lines
  for ii =1:1000
      # randomly get 4 candidate lines
      possibleK = 1:size(LINES,2);
      K = datasample(possibleK, 4,'Replace',false);
      lines = LINES(K);
      
      eq = zeros(3,length(lines));
      sumX = zeros(1, 4);
      sumY = zeros(1, 4);
      for k=1:4
          X1 = lines(k).point1'; #(x1,y1)'
          X2 = lines(k).point2'; #(x2,y2)'
          sumX(k) = X1(1) + X2(1);
          sumY(k) = X1(2) + X2(2);
      end
      
      # The top horizontal line
      i = find(sumY == min(sumY));
      i = i(randi(length(i)));
      X1 = lines(i).point1'; #(x1,y1)'
      X2 = lines(i).point2'; #(x2,y2)'
      eq(:,1) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
      
      # The right vertical line
      i = find(sumX == max(sumX));
      i = i(randi(length(i)));
      X1 = lines(i).point1'; #(x1,y1)'
      X2 = lines(i).point2'; #(x2,y2)'
      eq(:,2) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
      
      # The bottom horizontal line
      i = find(sumY == max(sumY));
      i = i(randi(length(i)));
      X1 = lines(i).point1'; #(x1,y1)'
      X2 = lines(i).point2'; #(x2,y2)'
      eq(:,3) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
      
      # The left vertical line
      i = find(sumX == min(sumX));
      i = i(randi(length(i)));
      X1 = lines(i).point1'; #(x1,y1)'
      X2 = lines(i).point2'; #(x2,y2)'
      eq(:,4) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
      
      
      # Find intersection of lines (corners) between the 4 random candidates
      points = zeros(2,4);   # first row: x-coordinate of the points, second row: y-coordinate
      for i = 1:1:length(eq)
          
          if i ~= length(eq)
              l1 = eq(:, i);
              l2 = eq(:, i+1);
              intersect = cross(l1, l2);
              intersect = intersect/intersect(3);
              points(:,i) = intersect(1:2);
          else
              l1 = eq(:, i);
              l2 = eq(:, 1);
              intersect = cross(l1, l2);
              intersect = intersect/intersect(3);
              points(:,i) = intersect(1:2);
          end
          
      end
      
      # check whether these points are better candidates for corners
      # area cannot be too big, all points must be within the image pixels
      area = polyarea(points(1,:), points(2,:));
      if area > max_area && ~isnan(area) && ~any(any((points<0))) && ...
              ~any(points(1,:)>size(image,2)) && ~any(points(2,:)>size(image,1)) 
          max_area = area;
          max_points = points;
      end
      
  end
  points = max_points;
      hold on
      for k = 1:length(LINES)
          xy = [LINES(k).point1; LINES(k).point2];
          plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
          # Plot beginnings and ends of lines
          plot(xy(1,1),xy(1,2),'x','LineWidth',2,'markersize', 20, 'Color','yellow');
          plot(xy(2,1),xy(2,2),'x','LineWidth',2,'markersize', 20, 'Color','red');
      end
      title('image with multiple line candidates')
      
      figure
      imshow(image)
      hold on
      plot(points(1,:),points(2,:),'rx','markersize', 30, 'linewidth',5)
      hold on
      plot(points(1,1:2),points(2,1:2),'Color','green','LineWidth',3)
      hold on
      plot(points(1,2:3),points(2,2:3),'Color','green','LineWidth',3)
      hold on
      plot(points(1,3:4),points(2,3:4),'Color','green','LineWidth',3)
      hold on
      plot([points(1,4),points(1,1)],[points(2,4),points(2,1)],...
          'Color','green','LineWidth',3)
      title('image with corners')
      '''
