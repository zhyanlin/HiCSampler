import cooler,sys
import numpy as np
cool = cooler.Cooler(sys.argv[1])


print(np.sqrt(np.nansum(cool.matrix(balance=False).fetch(sys.argv[2]))/np.nansum(cool.matrix(balance=True).fetch(sys.argv[2]))))
    