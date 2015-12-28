#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats


def qqplot(data,thedist,outfile,suptitle = ''):
	fig = plt.figure()
	fig.suptitle(suptitle,fontsize=14,fontweight='bold')
	stats.probplot(data,dist=thedist,plot=plt)	
	plt.savefig(outfile,format = 'pdf')

def qqplotpdf(data,thedist,outfile,suptitle = ''):
	outfilehandler = PdfPages(outfile)
	fig = plt.figure()
	fig.suptitle(suptitle,fontsize=14,fontweight='bold')
	stats.probplot(data,dist=thedist,plot=plt)	
	plt.savefig(outfilehandler,format = 'pdf')
	outfilehandler.close()
if __name__ == '__main__':
	outfile = PdfPages('qqplot.pdf')
	data = np.random.uniform(2,10,100)
	qqplot(data,'uniform',outfile)
	outfile.close()	

