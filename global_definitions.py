"""
@author: Michael Labenbacher
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys
paths = {'data':'data', 'plots':'plots','bg_removed':'bg_removed','bg_removed_RBM':'bg_removed_RBM'}
col_nm = {'633':'#ff4200',
          '568':'#ffe600',#'#c5e600',# 588 actually (looks better)
          '515':'#6EC610',
          '488':'#00c6cc',# 0.8 brightness adjusted
          '458':'#0071ff',
          '405':'#8200c8'}
col_rbm = {'1':'#000000','2':'#3A4785','3':'#C74336','4':'#6277DE'}
col_rbm_ = {'1':'k','2':'#943228','3':'#000000','4':'#3A4785','5':'#C74336'}
def col_adjust_brightness(color, amount=0.5): 
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

