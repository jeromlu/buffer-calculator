# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:37:42 2019

@author: jeromlu2
"""

from PyQt5.uic import loadUiType
 
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import numpy as np
from scipy import interpolate
import pandas as pd

    
Ui_MainWindow, QMainWindow = loadUiType('txfmr.ui')

class Main(QMainWindow, Ui_MainWindow):
    
    def __init__(self,):
        super(Main, self).__init__()
        self.setupUi(self)
        self.fstart = 2
        self.fstop = 4
        self.fstartLblLine.setText(str(self.fstart))
        self.fstopLblLine.setText(str(self.fstop))
        self.LpLblLine.setText('1')
        self.LsLblLine.setText('1')
        self.kLblLine.setText('0.75')
        self.w = np.linspace(2*np.pi*self.fstart*1e9,2*np.pi*self.fstop*1e9,101)
        self.DrawBtn.clicked.connect(self.draw)

    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas, self.mplwindow, coordinates=True)
        self.mplvl.addWidget(self.toolbar)
 
    def rmmpl(self,):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()
 
    def R(self,L,w,rf=1,rf2=1,deg=0):
        #rf = 1.78/1.24
        r0 = rf*L*1e9
        rf2 = 1
        alpha = rf2/(w[-1]-w[0])
        if deg == 0:
            return np.array([r0])
        return r0 + (alpha*(w-w[0]))**deg
    
    
    def calc_txf_std(self,Lp,Ls,k,w,rf=1,rf2=1,deg=0):
        RL = 50
        Rp = self.R(Lp,self.w,1.78/1.24,rf2,deg)                         #primary inductance resistance
        Rs = self.R(Ls,self.w,4.66/2.4,rf2,deg)
        tratio = 1/k*np.sqrt(Lp/Ls)
        n = k*np.sqrt(Lp/Ls)
        M = k*np.sqrt(Lp*Ls)
        Qip = w*Lp/Rp
        Qis = w*Ls/Rs
        Ql = w*Ls/(Rs+RL)
        Rsp = Rp + Ql**2/(1+Ql**2)*n**2*(Rs+RL)
        Lsp = Lp*(1-k**2*Ql**2/(1+Ql**2))
        Qsp = w*Lsp/Rsp
        Rpp = Rsp*(1+Qsp**2)
        Lpp = Lsp*(1+Qsp**2)/Qsp**2
        Cres = 1/(w**2*Lpp)
        Qc = 1/(w*Rp*Cres)
        w0 = 1/np.sqrt(Cres*Lp)
        ww0 = np.array([(w/wx)**2-1 for wx in w0])
        R_e_o = Rs + (w*M)**2/(Rp*(1+(Qc*ww0)**2))
        R_e_o_max = np.array([x.max() for x in R_e_o[:]])
        ww0 = (w/w0)**2-1
        IL = 10*np.log10(Rsp/(Ql**2/(1+Ql**2)*n**2*RL))
        Lprime_p = Lp*(1 - (1-k**2*Ql**2/(1+Ql**2))*(1+Qsp**2)/Qsp**2)
        Qprime_p = w*Lprime_p/Rp
        R_o_res = Rs + (w*M)**2/(Rp*(1+Qprime_p**2))
        X_o_res = w*Ls - (w*M)**2/(w*Lprime_p*(1+Qprime_p**2)/Qprime_p**2)
        RTL = 20*np.log10(abs((R_e_o_max-50)/(R_e_o_max+50)))
        return({'tratio':tratio,'n':n,'Ql':Ql,'Qsp':Qsp,'Rsp':Rsp,'Lsp':Lsp,'Rpp':Rpp,'Lpp':Lpp,'Cres':Cres,'IL':IL,\
               'Lprime_p':Lprime_p,'Qprime_p':Qprime_p,'R_o_res':R_o_res,'X_o_res':X_o_res,'RL':RTL,'Rs':Rs,'Rp':Rp,'Qip':Qip,'Qis':Qis})

    def plt_stuff(self,Lp=1,Ls=0.9,k=0.85):

        fig = Figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)

        ax1.plot(self.w*1e-9/(2*np.pi),self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['Rpp'],label='Rpp')
        ax1.legend()
        fig.suptitle('voltage ratio :{:.2f}'.format(self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['tratio']))
        ax1.grid(True)
        ax3.plot(self.w*1e-9/(2*np.pi),self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['IL'],label='IL (dB)',color='g')
        ax3.legend()
        ax3.grid(True)
        ax4.plot(self.w*1e-9/(2*np.pi),self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['RL'],label='RL (dB)',color='k')
        ax4.legend()
        ax4.grid(True)
        if len(self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['Rs']) == 1:
            ax2.plot([self.w[0]*1e-9/(2*np.pi),self.w[-1]*1e-9/(2*np.pi)],[self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['Rs'],self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['Rs']],label='Rs',color='r')
            ax2.plot([self.w[0]*1e-9/(2*np.pi),self.w[-1]*1e-9/(2*np.pi)],[self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['Rp'],self.calc_txf_std(Lp*1e-9,Ls*1e-9,k,self.w)['Rp']],label='Rp',color='b')
    
        else:
            ax2.plot(w*1e-9/(2*pi),calc_txf_std(Lp*1e-9,Ls*1e-9,k,w)['Rs'],label='Rs',color='r')
            ax2.plot(w*1e-9/(2*pi),calc_txf_std(Lp*1e-9,Ls*1e-9,k,w)['Rp'],label='Rp',color='b')
        ax2.legend()
        ax2.grid(True)
        self.rmmpl()
        self.addmpl(fig)

    def draw(self,):
        Lp = float(self.LpLblLine.text())
        Ls= float(self.LsLblLine.text())
        k= float(self.kLblLine.text())
        self.fstart = float(self.fstartLblLine.text())
        self.fstop = float(self.fstopLblLine.text())
        self.w = np.linspace(2*np.pi*self.fstart*1e9,2*np.pi*self.fstop*1e9,101)
        self.plt_stuff(Lp,Ls,k)


if __name__ == '__main__':
    import sys
    from PyQt5 import QtGui, QtWidgets
    import numpy as np
    import pandas as pd
    fig1 = Figure()
    app = QtWidgets.QApplication(sys.argv)
    main = Main()
    main.addmpl(fig1)
    main.show()
    sys.exit(app.exec_())


