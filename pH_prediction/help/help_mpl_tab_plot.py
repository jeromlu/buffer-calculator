# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:40:30 2019

@author: jeromlu2
"""
#third party packages
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QAction
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem, QVBoxLayout
from PyQt5.QtWidgets import QHeaderView, QPushButton, QLineEdit, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QTabWidget, QFrame
from PyQt5.QtWidgets import QAbstractScrollArea, QDialog


from PyQt5.QtCore import Qt


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure




import itertools

class MultiTabNavTool(NavigationToolbar):
    #====================================================================================================
    def __init__(self, canvases, tabs, parent=None):
        self.canvases = canvases
        self.tabs = tabs

        NavigationToolbar.__init__(self, canvases[0], parent)

    #====================================================================================================
    def get_canvas(self):
        return self.canvases[self.tabs.currentIndex()]

    def set_canvas(self, canvas):
        self._canvas = canvas

    canvas = property(get_canvas, set_canvas)
    
class MultiTabNavTool_v2(QWidget):

    def __init__(self, canvases, tabs, parent=None):
        QWidget.__init__(self, parent)
        self.canvases = canvases
        self.tabs = tabs
        self.toolbars = [NavigationToolbar(canvas, parent) for canvas in self.canvases]
        vbox = QVBoxLayout()
        for toolbar in self.toolbars:
            vbox.addWidget(toolbar)
        self.setLayout(vbox)
        self.switch_toolbar()
        self.tabs.currentChanged.connect(self.switch_toolbar)

    def switch_toolbar(self):
        for toolbar in self.toolbars:
            toolbar.setVisible(False)
        self.toolbars[self.tabs.currentIndex()].setVisible(True)

class MplMultiTab(QMainWindow):
    #====================================================================================================
    def __init__(self, parent=None, figures=None, labels=None):
        QMainWindow.__init__(self, parent)

        self.main_frame = QWidget()
        self.tabWidget = QTabWidget( self.main_frame )
        self.create_tabs( figures, labels )

        # Create the navigation toolbar, tied to the canvas
        self.mpl_toolbar = MultiTabNavTool_v2(self.canvases, self.tabWidget, self.main_frame)

        self.vbox = vbox = QVBoxLayout()
        vbox.addWidget(self.mpl_toolbar)
        vbox.addWidget(self.tabWidget)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    #====================================================================================================
    def create_tabs(self, figures, labels ):

        if labels is None:      labels = []
        figures =  [Figure()] if figures is None else figures     #initialise with empty figure in first tab if no figures provided
        self.canvases = [self.add_tab(fig, lbl) 
                            for (fig, lbl) in itertools.zip_longest(figures, labels) ]

    #====================================================================================================
    def add_tab(self, fig=None, name=None):
        '''dynamically add tabs with embedded matplotlib canvas with this function.'''

        # Create the mpl Figure and FigCanvas objects. 
        if fig is None:
            fig = Figure()
            ax = fig.add_subplot(111)

        canvas = fig.canvas if fig.canvas else FigureCanvas(fig)
        canvas.setParent(self.tabWidget)
        canvas.setFocusPolicy(Qt.ClickFocus)

        #self.tabs.append( tab )
        name = 'Tab %i'%(self.tabWidget.count()+1) if name is None else name
        self.tabWidget.addTab(canvas, name)

        return canvas
    
    def closeEvent(self, evt):
        QApplication.quit()
    
import numpy as np
import matplotlib.pyplot as plt
import sys

x = np.linspace(1, 2 * np.pi, 100)
figures = []
for i in range(1,5):
    fig, ax = plt.subplots()
    y = np.sin(np.pi*i*x)+0.1*np.random.randn(100)
    ax.plot(x,y)
    figures.append( fig )

app = QApplication(sys.argv)
ui = MplMultiTab( figures=figures )
ui.show()
app.exec_()