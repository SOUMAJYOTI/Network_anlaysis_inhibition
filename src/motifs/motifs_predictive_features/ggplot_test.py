from ggplot import *
# our trusty old friends
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

p = ggplot(aes(x='wt', y='mpg'), data=mtcars)
print(p + geom_point() + facet_grid("cyl", "gear", scales="free_y"))




