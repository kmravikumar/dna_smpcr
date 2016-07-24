#!/usr/bin/python

# plotting template and primer 
import matplotlib.pyplot as plt

def plot_template(template,primer=0,yt=1):
    """
    plot template using a number line
    if primer then primer flag should be 1
    set yt for where to plot the line 
    """
    plt_min = template[0] - 2
    plt_max = template[-1] + 2
    plt.xlim((plt_min,plt_max))
    plt.ylim((yt-1,2))
    if primer==1:
        plt_min = plt_min + 2 -0.2
    x = np.arange(plt_min,plt_max,1)
    y = yt*np.ones(len(x))
    plt.plot(x,y,'-k',lw=2)
    y2 = yt*np.ones(len(template))
    plt.plot(template,y2,'.k',ms=30)
    plt.axis('off')
    
