#!/usr/bin/env python
import numpy as np
import scipy
import matplotlib.pylab as plt
import sys
import scanner

def main(iname):
    image = plt.imread("./pictures/"+iname+".jpg")
    fig, axs = plt.subplots(1,4,sharex=True, sharey=True, figsize=(19.5,10))
    axs[0].imshow(image)
    scanner.scanner(image, True, axs)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    for a in sys.argv[1:]:
      main(a)
