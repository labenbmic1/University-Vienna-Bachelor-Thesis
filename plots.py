"""
@author: Michael Labenbacher
"""
import numpy as np
#import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
#import matplotlib.patches as pltp
plt.style.use("labor.mplstyle")
plt.rcParams["text.latex.preamble"] = r"""
\usepackage{amsmath}\usepackage{amssymb}\usepackage{siunitx}
"""
from matplotlib.ticker import FormatStrFormatter
#from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from matplotlib.legend import Legend
#from matplotlib.legend_handler import update_from_first_child, HandlerPatch
import convert_data
import data_preprocessing
#import model
# define plot colors
import global_definitions

def color_selector(data):
    col = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        wavelength, file_number = data_preprocessing.get_wavelength_from_label(data_key, file_number=True)
        quantity = data_preprocessing.get_filequantity_for_wavelength(data,wavelength)
        file_step_size = 0.2 if quantity <= 4 else (0.05 if quantity <= 10 else (0.025 if quantity < 20 else 0.01))
        col[i] = 'b'
        if wavelength in global_definitions.col_nm:
            col[i] = global_definitions.col_adjust_brightness(global_definitions.col_nm[wavelength], np.max([0,1-file_step_size*file_number]))
    return col

def rescale_and_clipping(x,y,quantity_datas,x_range=(0,5000),rescale='global'):
    y_max = 0
    
    if rescale == 'global':
        for i in range(quantity_datas):
            x_min_idx = np.argmax(x[i]>=x_range[0])
            x_max_idx = np.argmax(x[i]>=x_range[1])
            x[i] = x[i][x_min_idx:x_max_idx]
            y[i] = y[i][x_min_idx:x_max_idx]
            temp = np.max(y[i])
            if y_max < temp:
                y_max = temp
        for i in range(quantity_datas):
            y[i] = convert_data.max_to_one(y[i], y_max=y_max)
    else:
        for i in range(quantity_datas):
            x_min_idx = np.argmax(x[i]>=x_range[0])
            x_max_idx = np.argmax(x[i]>=x_range[1])
            if x[i][-1] <= x_range[1]:
                x_max_idx = len(x[i])-1
            x[i] = x[i][x_min_idx:x_max_idx]
            y[i] = y[i][x_min_idx:x_max_idx]
            y[i] = convert_data.max_to_one(y[i])
    return x,y,y_max

def rayleigh(res_fit=(np.zeros(1),np.zeros(1)), res=(np.zeros(1),np.zeros(1)), init_fit=(np.zeros(1),np.zeros(1)), init=(np.zeros(1),np.zeros(1)), plot_range=(0,2), label='633', **kwargs):
    # define figure
    fig, ax = plt.subplots(figsize=(6.472*1.25, 3.236*1.25))
    
    # data
    x_res,y_res = res
    x_res_fit,y_res_fit = res_fit
    x_init,y_init = init
    x_init_fit,y_init_fit = init_fit
    
    # color
    col = color_selector({
        label+'-1':np.column_stack((x_res, y_res)), label+'-2':np.column_stack((x_res_fit, y_res_fit)),
        label+'-3':np.column_stack((x_init, y_init)), label+'-4':np.column_stack((x_init_fit, y_init_fit))
        })
    
    # plot
    if 'y_max' not in kwargs:
        y_res = convert_data.max_to_one(y_res)
        y_res_fit = convert_data.max_to_one(y_res_fit)
        y_init = convert_data.max_to_one(y_init)
        y_init_fit = convert_data.max_to_one(y_init_fit)
    else:
        if kwargs['y_max'] is not None:
            y_res = convert_data.max_to_one(y_res, kwargs['y_max'])
            y_res_fit = convert_data.max_to_one(y_res_fit, kwargs['y_max'])
            y_init = convert_data.max_to_one(y_init, kwargs['y_max'])
            y_init_fit = convert_data.max_to_one(y_init_fit, kwargs['y_max'])

    ax.plot(x_res, y_res, marker='x', linestyle='', color=col[0], markersize=2)
    ax.plot(x_res_fit, y_res_fit, marker='', linestyle='-', color=col[1])
    ax.plot(x_init, y_init, marker='o', linestyle='', color=col[2], markersize=2)
    ax.plot(x_init_fit, y_init_fit, marker='', linestyle='-', color=col[3])
    
    plt.show()

def spectroscope_raw(data: dict, inset=[None], **kwargs):
    # define figure
    fig, ax = plt.subplots(figsize=(6.472*1.25, 3.236*1.25))
    
    if inset[0] is not None:
        temp_ = 1 if len(inset) == 1 else 2
        axin = inset_axes(ax, width=0.8*temp_, height=0.8*temp_, 
                          bbox_to_anchor=(0.05, 0, 1, 0.95),
                          bbox_transform=ax.transAxes, loc=2)
    
    # color
    col = color_selector(data)
    
    # plot
    i_inset = 0
    for i, (data_key, data_value) in enumerate(data.items()):
        # separate data and process
        x = data_value[:,0]
        y = data_value[:,1]
        if 'y_max' not in kwargs:
            if 'y_max_region' not in kwargs:
                y = convert_data.max_to_one(y)
            else:
                if kwargs['y_max_region'][0] is not None:       
                    x_min_idx = np.argmax(x>=kwargs['y_max_region'][i][0])
                    x_max_idx = np.argmax(x>=kwargs['y_max_region'][i][1])
                    if x_max_idx == 0:
                        x_max_idx = kwargs['y_max_region'][i][1]
                    y = convert_data.max_to_one(y,y_max_region=[x_min_idx,x_max_idx])
        else:
            if kwargs['y_max'] is not None:
                y = convert_data.max_to_one(y, kwargs['y_max'][i])
        if i not in inset:
            ax.plot(x, y, marker='.', linestyle='-', color=col[i])
        else: 
            axin.plot(x,y+i_inset,marker='.', linestyle='-', color=col[i])
            i_inset += 1
            
            label = data_key[:3]
            y_pos = 0.3+(i_inset-1)
            if i_inset == 3:
                y_pos += 0.1
            if i_inset == 4:
                y_pos -= 0.07
            l_x = 270
            axin.annotate(label, xy=(l_x,y_pos), xytext=(l_x,y_pos), color=col[i], ha='center', va='center')
            
        
    
    if 'fits' in kwargs:
        for fit in kwargs['fits']:
            if 'y_max' not in kwargs:
                if len(fit) > 2 and 'y_max' in fit[2]:
                    fit[1] = convert_data.max_to_one(fit[1], fit[2]['y_max'])        
                else:
                    fit[1] = convert_data.max_to_one(fit[1])
            else:
                if kwargs['y_max'] is not None:
                    fit[1] = convert_data.max_to_one(fit[1], kwargs['y_max'][i])
            ax.plot(fit[0], fit[1], marker='', linestyle='--', color='k')
    
    # dummy for label
    for i, (data_key, data_value) in enumerate(data.items()):
        if i<9:
            ax.plot([],[], linestyle="-", marker='.', linewidth=1, markersize=5, color=col[i], label=data_key)

    ax.legend()
    
    # plot settings
    #ax.set_xlim(-100, 500)
    if 'bg_removal_range' in kwargs:
        ax.set_ylim(-0.05, 1)
        xmin = np.min(kwargs['bg_removal_range'])
        xmax = np.max(kwargs['bg_removal_range'])
        ax.set_xlim(xmin, xmax)
        ax.set_aspect((xmax-xmin)/(1+0.05)*0.5)
    else:
        ax.set_ylim(0, 1)
        ax.set_xlim(-50, 3500)
        ax.set_aspect((3500+100)/1*0.5)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # inset 
    if inset[0] is not None:
        axin.set_xlim(0, 300)
        axin.set_ylim(0, len(inset))
        #axin.set_yticks([i for i in range(80, 140+1, 20)])
        #axin.set_xticks([i*0.5 for i in range(0, 2*eta+1, 1)])
        axin.tick_params(labelleft=False, labelbottom=False)

    # legend
    ax_legend = ax.legend(loc='upper right')
    for legobj in ax_legend.legend_handles:
        legobj.set_linewidth(1)

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'spectroscope_raw'+'.png')
    plt.show()

def G(data: dict, wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = 1400
        x_max = 1700
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots()
    
    # color
    col = color_selector(data)
            
    # separate data and process
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)

    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i], marker='x', linestyle='-', label=wavelengths[i], color=col[i])
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    # plot settings
    ax.set_ylim(0, 1.05)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect((x_max - x_min)/1.05)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'G'+'.png')
    plt.show()
    
def G_split(data: dict, wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = 1400
        x_max = 1700
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots()
    
    # color
    col = color_selector(data)
            
    # separate data and process
    quantity_data = len(data)
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
        
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i]+i, marker='x', linestyle='-', label=wavelengths[i], color=col[i])
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    # plot settings
    ax.set_ylim(-0.15, quantity_data+0.15)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect((x_max - x_min)/(quantity_data+0.15))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'G_split'+'.png')
    plt.show()    

def DD(data: dict, wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = 2600
        x_max = 2800
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots()
    
    # color
    col = color_selector(data)
            
    # separate data and process
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
        
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i], marker='x', linestyle='-', label=wavelengths[i], color=col[i])
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    # plot settings
    ax.set_ylim(0, 1.05)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect((x_max - x_min)/1.05)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'DD'+'.png')
    plt.show()

def DD_split(data: dict, wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = 2600
        x_max = 2800
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots()
    
    # color
    col = color_selector(data)
            
    # separate data and process
    quantity_data = len(data)
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
        
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i]+i, marker='x', linestyle='-', label=wavelengths[i], color=col[i])
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    # plot settings
    ax.set_ylim(-0.15, quantity_data+0.15)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect((x_max - x_min)/(quantity_data+0.15))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'DD_split'+'.png')
    plt.show()   

def G_split_plt(data: dict, labels=[None,None], wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = 1400
        x_max = 1700
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots(figsize=(6.472*2, 3.236*2))
    
    # color
    col = color_selector(data)
            
    # separate data and process
    quantity_data = len(data)
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
        
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    di=0.2
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i]+i*(1+di), marker='x', linestyle='-', label=wavelengths[i], color=col[i])
        ax.plot(x[i], np.full(len(x[i]),i*(1+di)), marker='', linestyle='--', label=wavelengths[i], color='k', linewidth=0.2)
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    dd = [0]*len(labels[0])
    if labels[0] is not None:
        for i,(l_x, l_y) in enumerate(zip(labels[0],labels[1])):
            y_pos = i*(1+di)+0.17
            l_y = l_y + r'\,eV'
            ax.annotate(l_y, xy=(l_x+dd[i],y_pos), xytext=(l_x+dd[i],y_pos), color=col[i], ha='center', va='center')
    
    
    # plot settings
    ax.tick_params(axis='y',which='both',labelleft=False,left=False,right=False)
    ax.set_ylim(-0.15, quantity_data*(1+di)+0.15)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect(2*(x_max - x_min)/(quantity_data+0.15))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'G_split_plt'+'.png')
    plt.show()    

def DD_split_plt(data: dict, labels=[None,None], wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = 2600
        x_max = 2800
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots(figsize=(6.472*2, 3.236*2))
    
    # color
    col = color_selector(data)
            
    # separate data and process
    quantity_data = len(data)
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
        
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    di=0.2
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i]+i*(1+di), marker='x', linestyle='-', label=wavelengths[i], color=col[i])
        ax.plot(x[i], np.full(len(x[i]),i*(1+di)), marker='', linestyle='--', label=wavelengths[i], color='k', linewidth=0.2)
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    dd = [0]*len(labels[0])
    if labels[0] is not None:
        for i,(l_x, l_y) in enumerate(zip(labels[0],labels[1])):
            y_pos = i*(1+di)+0.2 
            if i==0:
                y_pos += 0.1
            l_y = l_y + r'\,eV'
            ax.annotate(l_y, xy=(l_x+dd[i],y_pos), xytext=(l_x+dd[i],y_pos), color=col[i], ha='center', va='center')
    
    
    # plot settings
    ax.tick_params(axis='y',which='both',labelleft=False,left=False,right=False)
    ax.set_ylim(-0.15, quantity_data*(1+di)+0.15)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect(2*(x_max - x_min)/(quantity_data+0.15))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'DD_split_plt'+'.png')
    plt.show()   
    
    
def DD_split_plt_(data: dict, dummy=[None,None],labels=[None,None], wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = 2600
        x_max = 2800
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots(figsize=(6.472*2, 3.236*2))
    
    # color
    col = color_selector(data)
            
    # separate data and process
    quantity_data = len(data)
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
        
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    di=0.2
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i]+i*(1+di), marker='x', linestyle='-', label=wavelengths[i], color=col[i])
        ax.plot(x[i], np.full(len(x[i]),i*(1+di)), marker='', linestyle='--', label=wavelengths[i], color='k', linewidth=0.2)
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    dd = [0]*len(labels[0])
    if labels[0] is not None:
        for i,(l_x, l_y) in enumerate(zip(labels[0],labels[1])):
            y_pos = i*(1+di)+0.2 
            if i==0:
                y_pos += 0.1
            l_y = l_y + r'\,eV'
            ax.annotate(l_y, xy=(l_x+dd[i],y_pos), xytext=(l_x+dd[i],y_pos), color=col[i], ha='center', va='center')
    
    j=0
    for i,d in enumerate(dummy[1]):
        ax.plot(dummy[0], d+j*(1+di), marker='', linestyle='-.', color='k')
        if i%2!=0:
            j+=1
    
    # plot settings
    ax.tick_params(axis='y',which='both',labelleft=False,left=False,right=False)
    ax.set_ylim(-0.15, quantity_data*(1+di)+0.15)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect(2*(x_max - x_min)/(quantity_data+0.15))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'DD_split_plt_'+'.png')
    plt.show()   
    
def RBM_split_plt(data: dict, labels=[None,None], wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = -50
        x_max = 300
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots(figsize=(6.472*2, 3.236*2))
    
    # color
    col = color_selector(data)
            
    # separate data and process
    quantity_data = len(data)
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
    
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']

    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)

    # plot
    di=0.2
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i]+i*(1+di), marker='x', linestyle='-', label=wavelengths[i], color=col[i])
        ax.plot(x[i], np.full(len(x[i]),i*(1+di)), marker='', linestyle='--', label=wavelengths[i], color='k', linewidth=0.2)
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    dd = [0]*len(labels[0])
    dy = [0.25,0.23,0.37,0.18]
    if labels[0] is not None:
        for i,(l_x, l_y) in enumerate(zip(labels[0],labels[1])):
            y_pos = i*(1+di)+dy[i]
            l_y = l_y + r'\,eV'
            ax.annotate(l_y, xy=(l_x+dd[i],y_pos), xytext=(l_x+dd[i],y_pos), color=col[i], ha='center', va='center')
    
    # plot settings
    ax.tick_params(axis='y',which='both',labelleft=False,left=False,right=False)
    ax.set_ylim(-0.15, quantity_data*(1+di)+0.15)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect(2*(x_max - x_min)/(quantity_data+0.15))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'RBM_split_plt'+'.png')
    plt.show()    

def RBM_split(data: dict, wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = -50
        x_max = 300
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots()
    
    # color
    col = color_selector(data)
            
    # separate data and process
    quantity_data = len(data)
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)
    
    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']

    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)

    # plot
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i]+i, marker='x', linestyle='-', label=wavelengths[i], color=col[i])
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    # plot settings
    ax.set_ylim(-0.15, quantity_data+0.15)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect((x_max - x_min)/(quantity_data+0.15))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'RBM_split'+'.png')
    plt.show()    
    
def RBM(data: dict, wavelengths=np.array([633, 568], dtype=np.float64), fits=[[None,None],[None,None]], **kwargs):    
    # plot range
    if 'plot_range' not in kwargs:
        x_min = -50
        x_max = 300
    else:
        x_min = kwargs['plot_range'][0]
        x_max = kwargs['plot_range'][1]
    
    # define figure
    fig, ax = plt.subplots()
    
    # color
    col = color_selector(data)
            
    # separate data and process
    x = [None]*len(data)
    y = [None]*len(data)
    for i, (data_key, data_value) in enumerate(data.items()):
        x[i] = np.array([d[0] for d in data_value], dtype=np.float64)
        y[i] = np.array([d[1] for d in data_value], dtype=np.float64)

    # rescale
    rescale='global'
    if 'rescale' not in kwargs:
        rescale = 'global'
    else:
        rescale = kwargs['rescale']
    x,y,y_max = rescale_and_clipping(x,y,len(data),x_range=(x_min,x_max),rescale=rescale)
    
    # plot
    for i, (data_key, data_value) in enumerate(data.items()):
        ax.plot(x[i], y[i], marker='x', linestyle='-', label=wavelengths[i], color=col[i])
    if fits[0][0] is not None:
        print("Not implemented currently")
    
    # plot settings
    ax.set_ylim(0, 1.05)
    ax.set_xlim(x_min, x_max)
    ax.set_aspect((x_max - x_min)/1.05)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.set_xlabel(r'Raman frequency in \unit{\per\centi\metre}')
    ax.set_ylabel(r'Intensity in a.u.')

    # save 
    plt.savefig(global_definitions.paths['plots']+'/'+'G'+'.png')
    plt.show()