#CREATION OF JPDF DATABASE

import sys
from create_h_and_h_wc import *


#relative_mld = True : 
#gap_list --> depth = MLD + gap_list[i]*MLD
# -1 --> z = Surface
# 1 --> z = 2*MLD

#relative_mld = False : 
#gap_list = depth

month = sys.argv[1]

gap_list = np.arange(0,620,30) #tous les 30m jusqu'a 600 m 
tstart = 0
tend = 230

(H_jpdf,H_wc) = create_H_nc(month,0,tend,gap_list,relative_mld=False)






