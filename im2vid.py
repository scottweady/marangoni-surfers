import numpy as np
from natsort import natsorted
import cv2
import glob
import sys

print("loading files")
filepath = sys.argv[1]
videoname = sys.argv[2]
filenames = glob.glob(filepath + "/*.png")
filenames = natsorted(filenames)
N = len(filenames)
img = []

for n in range(0,N):
  img.append(cv2.imread(filenames[n]))

height,width,layers=img[0].shape
fourcc = cv2.VideoWriter_fourcc(*'MP4V')
video = cv2.VideoWriter(videoname + '.mp4',fourcc,24,(width,height))

for n in range(0,N):
  print("writing file " + str(n) + " of " + str(N),end="\r")
  video.write(img[n])

cv2.destroyAllWindows()
video.release()

