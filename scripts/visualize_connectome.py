import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

in_csv = snakemake.input['connectome']
out_png = snakemake.output[0]

#load csv
C = np.loadtxt(in_csv)

#symmetrize C
C = C + C.transpose()

#plot it as image and save to file
plt.imshow(np.log(C+0.0001))
plt.savefig(out_png)

