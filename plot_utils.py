import matplotlib.pyplot as plt
import numpy as np

def plot_vstats_anom(path,size,title,
					 time,z,
					 anom_data,
					 ylim,vmin, vmax,
					 cbar_label):
	fig, ax = plt.subplots(1,figsize=size)
	ax.set_title(title,fontsize=20)
	ax.set_ylim(ylim)
	im = ax.pcolormesh(time,z,anom_data.transpose())
	im.set_clim(vmin=vmin, vmax=vmax)
	cbar = fig.colorbar(im , ax = ax)
	cbar.ax.set_ylabel(cbar_label, size=20,rotation=270)
	cbar.ax.get_yaxis().labelpad = 25
	cbar.ax.tick_params(labelsize=15)
	ax.set_ylabel('Z [km]',size=20)
	ax.set_xlabel('Days since C0$_2$ change',size=20)
	ax.tick_params(axis='x',labelsize=15)
	ax.tick_params(axis='y',labelsize=15)
	return (fig,ax)
	
					 