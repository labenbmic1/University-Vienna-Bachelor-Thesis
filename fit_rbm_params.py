"""
@author: micha
"""

import numpy as np

inits=dict({'488-3':[dict({'a0':5e-2}),
        
        dict({'c':0.03,'x0':250,'gamma':5}),
        dict({'c':0.03,'x0':255,'gamma':5}),
        dict({'c':0.03,'x0':260,'gamma':5}),
        dict({'c':0.03,'x0':280,'gamma':5}),

        # 11 in 170 - 225
        dict({'c':0.04,'x0':165+5,'gamma':5}),
        dict({'c':0.04,'x0':169+2,'gamma':5}),
        dict({'c':0.04,'x0':171+1,'gamma':5}),
        dict({'c':0.04,'x0':175,'gamma':5}),
        dict({'c':0.04,'x0':178,'gamma':5}),
        dict({'c':0.04,'x0':180,'gamma':5}),
        dict({'c':0.04,'x0':200,'gamma':5}),
        dict({'c':0.04,'x0':202,'gamma':5}),
        dict({'c':0.04,'x0':204,'gamma':5}),
        dict({'c':0.04,'x0':216,'gamma':5}),
        dict({'c':0.04,'x0':222,'gamma':5}),

        # ~17 in 125 - 150 and inf below
        dict({'c':0.04,'x0':143,'gamma':5}),
        dict({'c':0.04,'x0':142,'gamma':5}),
        dict({'c':0.04,'x0':141,'gamma':5}),
        dict({'c':0.04,'x0':140,'gamma':5}),
        dict({'c':0.04,'x0':134,'gamma':5}),
        dict({'c':0.04,'x0':132,'gamma':5}),
        dict({'c':0.04,'x0':128,'gamma':5}),
        dict({'c':0.04,'x0':82,'gamma':5}),
        dict({'c':0.04,'x0':85,'gamma':5}),
        dict({'c':0.04,'x0':95,'gamma':5}),
        dict({'c':0.04,'x0':100,'gamma':5}),
        dict({'c':0.04,'x0':102,'gamma':5}),
        dict({'c':0.04,'x0':104,'gamma':5}),
        dict({'c':0.04,'x0':106,'gamma':5}),
        dict({'c':0.04,'x0':108,'gamma':5}),
        dict({'c':0.04,'x0':110,'gamma':5}),
        dict({'c':0.04,'x0':112,'gamma':5})
    ],
    '633-3':[
        dict({'a0':5e-2}),
        #dict({'c':0.3,'x0':90,'sigma':10}),
        
        dict({'c':0.03,'x0':250,'gamma':5}),
        dict({'c':0.03,'x0':255,'gamma':5}),
        dict({'c':0.03,'x0':260,'gamma':5}),
        #dict({'c':0.03,'x0':270,'gamma':5}),
        dict({'c':0.03,'x0':280,'gamma':5}),
        # 17 in 170 - 225
        dict({'c':0.04,'x0':165+5,'gamma':5}),
        dict({'c':0.04,'x0':169+2,'gamma':5}),
        dict({'c':0.04,'x0':171+1,'gamma':5}),
        dict({'c':0.04,'x0':175,'gamma':5}),
        
        dict({'c':0.04,'x0':178,'gamma':5}),
        dict({'c':0.04,'x0':180,'gamma':5}),
        #dict({'c':0.04,'x0':182,'gamma':5}),
        #dict({'c':0.04,'x0':185,'gamma':5}),
        
        #dict({'c':0.04,'x0':190,'gamma':5}),
        #dict({'c':0.04,'x0':194,'gamma':5}),
        #dict({'c':0.04,'x0':197,'gamma':5}),
        dict({'c':0.04,'x0':200,'gamma':5}),
        dict({'c':0.04,'x0':202,'gamma':5}),
        dict({'c':0.04,'x0':204,'gamma':5}),
        #dict({'c':0.04,'x0':212,'gamma':5}),
        dict({'c':0.04,'x0':216,'gamma':5}),
        dict({'c':0.04,'x0':222,'gamma':5}),

        # ~14 in 125 - 150 and inf below
        #dict({'c':0.04,'x0':150,'gamma':5}),
        #dict({'c':0.04,'x0':148-4,'gamma':5}),
        dict({'c':0.04,'x0':143,'gamma':5}),
        dict({'c':0.04,'x0':142,'gamma':5}),
        dict({'c':0.04,'x0':141,'gamma':5}),
        dict({'c':0.04,'x0':140,'gamma':5}),
        #dict({'c':0.04,'x0':138,'gamma':5}),
        #dict({'c':0.04,'x0':136,'gamma':5}),
        dict({'c':0.04,'x0':134,'gamma':5}),
        dict({'c':0.04,'x0':132,'gamma':5}),
        #dict({'c':0.04,'x0':130,'gamma':5}),
        dict({'c':0.04,'x0':128,'gamma':5}),
        dict({'c':0.04,'x0':82,'gamma':5}),
        dict({'c':0.04,'x0':85,'gamma':5}),
        #dict({'c':0.04,'x0':88,'gamma':5}),
        #dict({'c':0.04,'x0':92,'gamma':5}),
        dict({'c':0.04,'x0':95,'gamma':5}),
        #dict({'c':0.04,'x0':100,'gamma':5}),
        #dict({'c':0.04,'x0':100,'gamma':5}),
        dict({'c':0.04,'x0':100,'gamma':5}),
        dict({'c':0.04,'x0':102,'gamma':5}),
        dict({'c':0.04,'x0':104,'gamma':5}),
        dict({'c':0.04,'x0':106,'gamma':5}),
        dict({'c':0.04,'x0':108,'gamma':5}),
        dict({'c':0.04,'x0':110,'gamma':5}),
        dict({'c':0.04,'x0':112,'gamma':5}),
    ]})
