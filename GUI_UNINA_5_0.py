"""
Created on Thu Apr 18 16:09:00 2019

@author: Alessandro Cianciulli, Il_Ciancio
@email: al.cianciulli@studenti.unina.it
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (QApplication, QWidget, QToolTip, QPushButton, QMessageBox)
from PyQt5.QtCore import QCoreApplication, Qt, QEvent
import numpy as np
import pyqtgraph as pg
from function import getDataForGUI, getDataFromChannel, getDataFromDatabase, Bmagnmirnv, tomo_centroid
from numpy import pi
import sys
import os
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib as mpl
import ipfn_rc
import dieti_rc
from scipy import linalg


def select_channel(name_data):
    if name_data == "Primary Current":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_093'
    elif name_data == "Horizontal Current":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_091'
    elif name_data == "Vertical Current":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_092'
    elif name_data == "Plasma Current":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_228'
    elif name_data == "Rc_Kalman":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_243'
    elif name_data == "Zc_Kalman":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_244'
    elif name_data == "Rc_Mirnov":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_083'
    elif name_data == "Zc_Mirnov":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_084'
    elif name_data == "Rc_Probes":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_081'
    elif name_data == "Zc_Probes":
        channel = 'MARTE_NODE_IVO3.DataCollection.Channel_082'
    else:
        channel = "NameError"
    return channel

class MyMplCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        self.compute_initial_figure()
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


# not avaible
class DynamicTomographyCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    global x_vessel, y_vessel
    N = 100
    th_vessel = np.linspace(0, 2*pi, N)
    x_vessel = 0.09 * np.cos(th_vessel)
    y_vessel = 0.09 * np.sin(th_vessel)

    def __init__(self, *args, **kwargs):
        global i1, x1, y1, t1
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)
        i1 = 0
        x1 = np.load("R_fromTomography.npy")
        y1 = np.load("z_fromTomography.npy")
        t1 = np.load("time_out.npy")

    def compute_initial_figure(self):
        global x_vessel, y_vessel
        self.axes.plot(x_vessel, y_vessel, color='k', linewidth=2.0)
        self.axes.set_title("VACUUM VESSEL")

    def update_figure(self):
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        global i1, x1, y1, t1, x_vessel, y_vessel
        self.axes.cla()
        self.axes.set_xlim([-0.1, 0.1])
        self.axes.set_ylim([-0.1, 0.1])
        if x1[i1] >= 0.085 or y1[i1] >= 0.085 or x1[i1] <= -0.085 or y1[i1] <= -0.085:
            self.axes.text(-0.02, 0, "NO PLASMA", bbox=dict(facecolor='red', alpha=0.5, pad=10))    
        self.axes.text(-0.1, 0.1, "time: "+str(t1[i1])+" microS")
        self.axes.plot(x1[i1], y1[i1], '-o', color='k')
        self.axes.plot(x_vessel, y_vessel, color='k', linewidth=2.0)
        self.axes.set_title("VACUUM VESSEL")
        self = self.draw()
        i1 = i1 + 1


class DynamicMirnovCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        global i2, x2, y2, t2
        i2 = 0
        x2 = np.load("R_fromMirnov.npy")
        y2 = np.load("z_fromMirnov.npy")
        t2 = np.load("time_out.npy")
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)

    def compute_initial_figure(self):
        global x_vessel, y_vessel
        self.axes.plot(x_vessel, y_vessel, color='b', linewidth=2.0)

    def update_figure(self):
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        global i2, x2, y2, t2, x_vessel, y_vessel
        self.axes.cla()
        self.axes.set_xlim([-0.1, 0.1])
        self.axes.set_ylim([-0.1, 0.1])
        if x2[i2] >= 0.085 or y2[i2] >= 0.085 or x2[i2] <= -0.085 or y2[i2] <= -0.085:
            self.axes.text(-0.02, 0, "NO PLASMA", bbox=dict(facecolor='red', alpha=0.5, pad=10))    
        self.axes.text(-0.1, 0.1, "time: "+str(t2[i2])+" microS")
        self.axes.plot(x2[i2], y2[i2], '-o', color='b')
        self.axes.plot(x_vessel, y_vessel, color='b', linewidth=2.0)
        self.axes.set_title("VACUUM VESSEL")
        self = self.draw()
        i2 = i2 + 1


class DynamicProbesCanvas(MyMplCanvas):
    def __init__(self, *args, **kwargs):
        global i3, x3, y3, t3
        i3 = 0
        x3 = np.load("R_fromProbes.npy")
        y3 = np.load("z_fromProbes.npy")
        t3 = np.load("time_out.npy")
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)

    def compute_initial_figure(self):
        global x_vessel, y_vessel
        self.axes.plot(x_vessel, y_vessel, color='r', linewidth=2.0)

    def update_figure(self):
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        global i3, x3, y3, t3, x_vessel, y_vessel
        self.axes.cla()
        self.axes.set_xlim([-0.1, 0.1])
        self.axes.set_ylim([-0.1, 0.1])
        if x3[i3] >= 0.085 or y3[i3] >= 0.085 or x3[i3] <= -0.085 or y3[i3] <= -0.085:
            self.axes.text(-0.02, 0, "NO PLASMA", bbox=dict(facecolor='red', alpha=0.5, pad=10))    
        self.axes.text(-0.1, 0.1, "time: "+str(t3[i3])+" microS")
        self.axes.plot(x3[i3], y3[i3], '-o', color='r')
        self.axes.plot(x_vessel, y_vessel, color='r', linewidth=2.0)
        self.axes.set_title("VACUUM VESSEL")
        self = self.draw()
        i3 = i3 + 1

class DynamicMinovVsPostionCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    global x1_vessel, y1_vessel
    N = 100
    th_vessel = np.linspace(0, 2*pi, N)
    x1_vessel = 9 * np.cos(th_vessel) + 46 
    y1_vessel = 9 * np.sin(th_vessel)

    def __init__(self, *args, **kwargs):
        global xT, yT, i4, xR, yR, t4, xP, yP, B_Mirnov, R_filaments, z_filaments, R_mirnov, z_mirnov, I_filament_all, Ip_magn_corr_value
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)
        N_FILAMENTS = 12
        i4 = 0
        xR = np.load("R_fromMirnov.npy")
        yR = np.load("z_fromMirnov.npy")
        xP = np.load("R_fromProbes.npy")
        yP = np.load("z_fromProbes.npy")
        xT = np.load("R_fromTomography.npy")
        yT = np.load("z_fromTomography.npy")
        t4 = np.load("time_out.npy")
        B_Mirnov = np.load("B_mirnov.npy")
        Ip_magn_corr_value = np.load("I_plasma.npy")
        xP = (1e02 * xP)+ 46 # in cm
        yP = (1e02 * yP)# in cm
        xR = (1e02 * xR)+ 46# in cm
        yR = (1e02 * yR)# in cm
        xT = (1e02 * xT)+ 46# in cm
        yT = (1e02 * yT)# in cm        
        #MIRNOV POSITION IN DEGREE
        R_mirnov = np.array([], np.float32)
        z_mirnov = np.array([], np.float32)
        ang_mirnov = -15
        for i in range(12):
            R_mirnov = np.append( R_mirnov, 9.35 * np.cos( np.deg2rad(ang_mirnov)  ) + 46 ) 
            z_mirnov = np.append( z_mirnov, 9.35 * np.sin( np.deg2rad(ang_mirnov) ) )
            ang_mirnov = ang_mirnov - 30
        # Lets draw the plasma filaments	
        R_filaments = np.array([], np.float32)
        z_filaments = np.array([], np.float32)
        degr_filament = 0
        radius = 5.5  # in [cm] (distance from the center of the chamber to the filaments)
        degr_fact = 360/N_FILAMENTS
        for i in range(0,N_FILAMENTS) :
            R_filaments = np.append( R_filaments,(46) + radius * np.cos( np.deg2rad(degr_filament) ))
            z_filaments = np.append( z_filaments, radius * np.sin( np.deg2rad(degr_filament) ) )
            degr_filament = degr_filament + degr_fact;
        # Matrix whose elements gives the contribution  to the measuremnt i  to a unitary current in the filament j [T]
        Mfp = np.zeros((12,N_FILAMENTS), np.float32)
        for i in range(12):
            for j in range(N_FILAMENTS):
                Mfp[i][j] = Bmagnmirnv(z_filaments[j], R_filaments[j], 1, R_mirnov[i], z_mirnov[i])

        Mpf = np.zeros((12,N_FILAMENTS), np.float32)
        Mpf = linalg.pinv(Mfp)
        I_filament_all = np.dot(Mpf,B_Mirnov)

    def compute_initial_figure(self):
        global xT, yT, i4, xR, yR, t4, xP, yP, B_Mirnov, R_filaments, z_filaments, R_mirnov, z_mirnov, I_filament_all, Ip_magn_corr_value
        self.axes.plot(x1_vessel, y1_vessel, color='k', linewidth=2.0)
        self.axes.set_title("VACUUM VESSEL")

    def update_figure(self):
        global xT, yT, i4, xR, yR, t4, xP, yP, B_Mirnov, R_filaments, z_filaments, R_mirnov, z_mirnov, I_filament_all, Ip_magn_corr_value
        self.axes.cla()
        name_current = 'Plasma Current: '+str(Ip_magn_corr_value[i4]) + ' [A]'
        name_time = 'Time: '+str(t4[i4])+' [ms]'
        self.axes.text(38,12,name_current, horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='red', alpha=0.5))
        self.axes.text(52,12,name_time, horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='red', alpha=0.5))
        #colors = I_filament_all[:,i4]
        abs_mirnov = np.abs(B_Mirnov)
        colors_mirnov = abs_mirnov[:,i4]
#        sc1 = self.axes.scatter(R_filaments, z_filaments, s=100, c = colors, cmap='autumn_r')
#        color = mpl.pyplot.colorbar(sc1, ax=self.axes)
#        color.set_label('Filaments Currrent [A]', rotation=270)
        sc2 = self.axes.scatter(R_mirnov, z_mirnov, marker='s',s=100, c = colors_mirnov, cmap='autumn_r',label='|Mirnov[T]|')
        color2 = mpl.pyplot.colorbar(sc2, ax=self.axes)
        color2.set_label('|Mirnov [T]|', rotation=270)
        for j in range(0, 12):
            self.axes.text(R_filaments[j],z_filaments[j]+1,str(int(I_filament_all[j,i4]))+" [A]", horizontalalignment='center', verticalalignment='center',bbox=dict(facecolor='red', alpha=0.5), label='Multifilament approssimation')
            self.axes.text(R_mirnov[j], z_mirnov[j],str(j+1))
            self.axes.plot(R_filaments[j],z_filaments[j], '.m', markersize=10)
        self.axes.plot(xP[i4], yP[i4], '.b', label='Centroid from probes' )
        self.axes.plot(xR[i4], yR[i4], '.r', label='Centroid from mirnov')
        self.axes.plot(xT[i4], yT[i4], '.k', label='Centroid from tomography')
        self.axes.plot(x1_vessel,y1_vessel,color='k', linewidth=2.0)
        self.axes.grid()
        self.axes.axis('equal')
        self.axes.set_xlabel('R[cm]')
        self.axes.set_ylabel('Z[cm]')
        self.axes.legend(loc='upper left')
        self = self.draw()
        #color.remove()
        color2.remove()
        i4 = i4 + 1

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1269, 886)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout.addWidget(self.label_2)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.gridLayout.addLayout(self.horizontalLayout, 0, 0, 2, 2)
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setMinimumSize(QtCore.QSize(156, 0))
        self.label_3.setMaximumSize(QtCore.QSize(156, 16777215))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignCenter)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 1, 2, 1, 2)
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setMovable(True)
        self.tabWidget.setObjectName("tabWidget")
        self.Centroid_position = QtWidgets.QWidget()
        self.Centroid_position.setObjectName("Centroid_position")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.Centroid_position)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_4 = QtWidgets.QLabel(self.Centroid_position)
        self.label_4.setMinimumSize(QtCore.QSize(152, 122))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.label_4.setFont(font)
        self.label_4.setAlignment(QtCore.Qt.AlignCenter)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)
        self.label_5 = QtWidgets.QLabel(self.Centroid_position)
        self.label_5.setMinimumSize(QtCore.QSize(152, 123))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.label_5.setFont(font)
        self.label_5.setAlignment(QtCore.Qt.AlignCenter)
        self.label_5.setObjectName("label_5")
        self.verticalLayout.addWidget(self.label_5)
        self.label_6 = QtWidgets.QLabel(self.Centroid_position)
        self.label_6.setMinimumSize(QtCore.QSize(152, 122))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.label_6.setFont(font)
        self.label_6.setAlignment(QtCore.Qt.AlignCenter)
        self.label_6.setObjectName("label_6")
        self.verticalLayout.addWidget(self.label_6)
        self.gridLayout_2.addLayout(self.verticalLayout, 0, 0, 1, 1)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        #PLOT CENTROID POSITION
        self.R_mirnov = pg.PlotWidget(self.Centroid_position) # QtWidgets.QGraphicsView(self.Centroid_position)
        self.R_mirnov.setMinimumSize(QtCore.QSize(639, 59))
        self.R_mirnov.setObjectName("R_mirnov")
        self.verticalLayout_2.addWidget(self.R_mirnov)
        self.Z_mirnov = pg.PlotWidget(self.Centroid_position) # QtWidgets.QGraphicsView(self.Centroid_position)
        self.Z_mirnov.setMinimumSize(QtCore.QSize(639, 58))
        self.Z_mirnov.setObjectName("Z_mirnov")
        self.verticalLayout_2.addWidget(self.Z_mirnov)
        self.R_probes = pg.PlotWidget(self.Centroid_position) # QtWidgets.QGraphicsView(self.Centroid_position)
        self.R_probes.setMinimumSize(QtCore.QSize(639, 58))
        self.R_probes.setObjectName("R_probes")
        self.verticalLayout_2.addWidget(self.R_probes)
        self.Z_probes =  pg.PlotWidget(self.Centroid_position) # QtWidgets.QGraphicsView(self.Centroid_position)
        self.Z_probes.setMinimumSize(QtCore.QSize(639, 58))
        self.Z_probes.setObjectName("Z_probes")
        self.verticalLayout_2.addWidget(self.Z_probes)
        self.R_tomgraphy = pg.PlotWidget(self.Centroid_position) # QtWidgets.QGraphicsView(self.Centroid_position)
        self.R_tomgraphy.setMinimumSize(QtCore.QSize(639, 58))
        self.R_tomgraphy.setObjectName("R_tomgraphy")
        self.verticalLayout_2.addWidget(self.R_tomgraphy)
        self.Z_tomography = pg.PlotWidget(self.Centroid_position) # QtWidgets.QGraphicsView(self.Centroid_position)
        self.Z_tomography.setMinimumSize(QtCore.QSize(639, 58))
        self.Z_tomography.setObjectName("Z_tomography")
        #END
        self.verticalLayout_2.addWidget(self.Z_tomography)
        self.gridLayout_2.addLayout(self.verticalLayout_2, 0, 1, 1, 1)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        #RZ ANIMATED PLOT
#        self.RZ_mirnov = QtWidgets.QGraphicsView(self.Centroid_position)
#        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
#        sizePolicy.setHorizontalStretch(0)
#        sizePolicy.setVerticalStretch(0)
#        sizePolicy.setHeightForWidth(self.RZ_mirnov.sizePolicy().hasHeightForWidth())
#        self.RZ_mirnov.setSizePolicy(sizePolicy)
#        self.RZ_mirnov.setObjectName("RZ_mirnov")
#        self.verticalLayout_3.addWidget(self.RZ_mirnov)
#        self.RZ_probes = QtWidgets.QGraphicsView(self.Centroid_position)
#        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
#        sizePolicy.setHorizontalStretch(0)
#        sizePolicy.setVerticalStretch(0)
#        sizePolicy.setHeightForWidth(self.RZ_probes.sizePolicy().hasHeightForWidth())
#        self.RZ_probes.setSizePolicy(sizePolicy)
#        self.RZ_probes.setObjectName("RZ_probes")
#        self.verticalLayout_3.addWidget(self.RZ_probes)
#        self.RZ_tomography = QtWidgets.QGraphicsView(self.Centroid_position)
#        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
#        sizePolicy.setHorizontalStretch(0)
#        sizePolicy.setVerticalStretch(0)
#        sizePolicy.setHeightForWidth(self.RZ_tomography.sizePolicy().hasHeightForWidth())
#        self.RZ_tomography.setSizePolicy(sizePolicy)
#        self.RZ_tomography.setObjectName("RZ_tomography")
#        self.verticalLayout_3.addWidget(self.RZ_tomography)
        #END
        self.gridLayout_2.addLayout(self.verticalLayout_3, 0, 2, 1, 1)
        self.tabWidget.addTab(self.Centroid_position, "")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_4.setObjectName("gridLayout_4")
        #MIRNOV VS CENTROID POSITION
#        self.graphicsView_5 = QtWidgets.QGraphicsView(self.tab)
#        self.graphicsView_5.setObjectName("graphicsView_5")
#        self.gridLayout_4.addWidget(self.graphicsView_5, 0, 0, 1, 1)
        #END
        self.tabWidget.addTab(self.tab, "")
        self.Database = QtWidgets.QWidget()
        self.Database.setObjectName("Database")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.Database)
        self.gridLayout_3.setObjectName("gridLayout_3")
        #PLOT FROM DATABASE
        self.graphicsView_3 = pg.PlotWidget(self.Database) # QtWidgets.QGraphicsView(self.Database)
        self.graphicsView_3.setObjectName("graphicsView_3")
        self.gridLayout_3.addWidget(self.graphicsView_3, 3, 0, 1, 1)
        self.graphicsView_4 =  pg.PlotWidget(self.Database) # QtWidgets.QGraphicsView(self.Database)
        self.graphicsView_4.setObjectName("graphicsView_4")
        self.gridLayout_3.addWidget(self.graphicsView_4, 4, 0, 1, 1)
        self.graphicsView =  pg.PlotWidget(self.Database) # QtWidgets.QGraphicsView(self.Database)
        self.graphicsView.setObjectName("graphicsView")
        self.gridLayout_3.addWidget(self.graphicsView, 1, 0, 1, 1)
        self.graphicsView_2 =  pg.PlotWidget(self.Database) # QtWidgets.QGraphicsView(self.Database)
        self.graphicsView_2.setObjectName("graphicsView_2")
        self.gridLayout_3.addWidget(self.graphicsView_2, 2, 0, 1, 1)
        #END
        self.Clear_button = QtWidgets.QPushButton(self.Database)
        self.Clear_button.setObjectName("Clear_button")
        self.gridLayout_3.addWidget(self.Clear_button, 0, 2, 1, 1)
        self.Plot_Button = QtWidgets.QPushButton(self.Database)
        self.Plot_Button.setObjectName("Plot_Button")
        self.gridLayout_3.addWidget(self.Plot_Button, 0, 1, 1, 1)   
        self.SelectChannel_combo = QtWidgets.QComboBox(self.Database)
        self.SelectChannel_combo.setObjectName("SelectChannel_combo")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.SelectChannel_combo.addItem("")
        self.gridLayout_3.addWidget(self.SelectChannel_combo, 1, 1, 1, 2)
        self.SelectChannel_combo_2 = QtWidgets.QComboBox(self.Database)
        self.SelectChannel_combo_2.setObjectName("SelectChannel_combo_2")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.SelectChannel_combo_2.addItem("")
        self.gridLayout_3.addWidget(self.SelectChannel_combo_2, 2, 1, 1, 2)
        self.SelectChannel_combo_3 = QtWidgets.QComboBox(self.Database)
        self.SelectChannel_combo_3.setObjectName("SelectChannel_combo_3")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.SelectChannel_combo_3.addItem("")
        self.gridLayout_3.addWidget(self.SelectChannel_combo_3, 3, 1, 1, 2)
        self.SelectChannel_combo_4 = QtWidgets.QComboBox(self.Database)
        self.SelectChannel_combo_4.setObjectName("SelectChannel_combo_4")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.SelectChannel_combo_4.addItem("")
        self.gridLayout_3.addWidget(self.SelectChannel_combo_4, 4, 1, 1, 2)
        self.tabWidget.addTab(self.Database, "")
        #new
                
        self.compare_tab = QtWidgets.QWidget()
        self.compare_tab.setObjectName("compare_tab")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.compare_tab)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.gridLayout_5 = QtWidgets.QGridLayout()
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.graphicsView_8 = pg.PlotWidget(self.compare_tab) # QtWidgets.QGraphicsView(self.compare_tab)
        self.graphicsView_8.setObjectName("graphicsView_8")
        self.gridLayout_5.addWidget(self.graphicsView_8, 3, 0, 1, 1)
        self.splitter = QtWidgets.QSplitter(self.compare_tab)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.widget = QtWidgets.QWidget(self.splitter)
        self.widget.setObjectName("widget")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.Plot_Button_Compare = QtWidgets.QPushButton(self.widget)
        self.Plot_Button_Compare.setObjectName("Plot_Button_Compare")
        self.verticalLayout_4.addWidget(self.Plot_Button_Compare)
        self.comboBox = QtWidgets.QComboBox(self.widget)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.verticalLayout_4.addWidget(self.comboBox)
        self.graphicsView_6 = pg.PlotWidget(self.compare_tab) # QtWidgets.QGraphicsView(self.widget)
        self.graphicsView_6.setObjectName("graphicsView_6")
        self.verticalLayout_4.addWidget(self.graphicsView_6)
        self.widget1 = QtWidgets.QWidget(self.splitter)
        self.widget1.setObjectName("widget1")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.widget1)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.Clear_button_Compare = QtWidgets.QPushButton(self.widget1)
        self.Clear_button_Compare.setObjectName("Clear_button_Compare")
        self.verticalLayout_5.addWidget(self.Clear_button_Compare)
        self.comboBox_6 = QtWidgets.QComboBox(self.widget1)
        self.comboBox_6.setObjectName("comboBox_6")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.verticalLayout_5.addWidget(self.comboBox_6)
        self.graphicsView_7 = pg.PlotWidget(self.compare_tab) # QtWidgets.QGraphicsView(self.widget1)
        self.graphicsView_7.setObjectName("graphicsView_7")
        self.verticalLayout_5.addWidget(self.graphicsView_7)
        self.gridLayout_5.addWidget(self.splitter, 2, 0, 1, 1)
        self.gridLayout_6.addLayout(self.gridLayout_5, 0, 0, 1, 1)
        self.tabWidget.addTab(self.compare_tab, "")
        
        #end
        self.gridLayout.addWidget(self.tabWidget, 3, 0, 6, 2)
        self.Start_Button = QtWidgets.QPushButton(self.centralwidget)
        self.Start_Button.setMinimumSize(QtCore.QSize(75, 0))
        self.Start_Button.setMaximumSize(QtCore.QSize(75, 16777215))
        self.Start_Button.setObjectName("Start_Button")
        self.gridLayout.addWidget(self.Start_Button, 3, 2, 1, 1)
        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setMinimumSize(QtCore.QSize(156, 0))
        self.textBrowser.setMaximumSize(QtCore.QSize(156, 16777215))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.textBrowser.setFont(font)
        self.textBrowser.setObjectName("textBrowser")
        self.textBrowser.append("\nSelect or Write a Shot Number and after press Start!")
        self.textBrowser.append("\nMAX_EVENT load the last shot available!")
        self.textBrowser.append("\nNegative value load the (MAX_EVENT-VALUE) shot!")
        self.textBrowser.append("\nTOMOGRAPHY is valid only for shots after #44880!") 
        self.textBrowser.append("\KALMAN is valid only for shots after #46881!")
        self.gridLayout.addWidget(self.textBrowser, 4, 2, 5, 2)
        self.Restart_button = QtWidgets.QPushButton(self.centralwidget)
        self.Restart_button.setMinimumSize(QtCore.QSize(75, 0))
        self.Restart_button.setMaximumSize(QtCore.QSize(75, 16777215))
        self.Restart_button.setObjectName("Restart_button")
        self.gridLayout.addWidget(self.Restart_button, 3, 3, 1, 1)
        self.ShotNumber_combo = QtWidgets.QComboBox(self.centralwidget)
        self.ShotNumber_combo.setMinimumSize(QtCore.QSize(156, 0))
        self.ShotNumber_combo.setMaximumSize(QtCore.QSize(156, 16777215))
        self.ShotNumber_combo.setEditable(True)
        self.ShotNumber_combo.setObjectName("ShotNumber_combo")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.ShotNumber_combo.addItem("")
        self.gridLayout.addWidget(self.ShotNumber_combo, 2, 2, 1, 2)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "CIANCIO'S GUI"))
        self.label_2.setText(_translate("MainWindow", "<html><head/><body><p><img src=\":/photo_ipfn/IPFN_color.png\"/></p></body></html>"))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p><img src=\":/dieti/DIETI.png\"/></p></body></html>"))
        self.label_3.setText(_translate("MainWindow", "SHOT NUMBER"))
        self.label_4.setText(_translate("MainWindow", "MIRNOV COILS"))
        self.label_5.setText(_translate("MainWindow", "ELECTRIC PROBES"))
        self.label_6.setText(_translate("MainWindow", "TOMOGRAPHY"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Centroid_position), _translate("MainWindow", "Centroid Position"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Mirnov Coils vs Centroid Position"))
        self.Clear_button.setText(_translate("MainWindow", "CLEAR"))
        self.Plot_Button.setText(_translate("MainWindow", "PLOT"))
        self.Plot_Button_Compare.setText(_translate("MainWindow", "PLOT"))
        self.Clear_button_Compare.setText(_translate("MainWindow", "CLEAR"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.compare_tab), _translate("MainWindow", "Compare"))
        self.SelectChannel_combo.setItemText(0, _translate("MainWindow", "Plasma Current"))
        self.SelectChannel_combo.setItemText(1, _translate("MainWindow", "Primary Current"))
        self.SelectChannel_combo.setItemText(2, _translate("MainWindow", "Horizontal Current"))
        self.SelectChannel_combo.setItemText(3, _translate("MainWindow", "Vertical Current"))
        self.SelectChannel_combo.setItemText(4, _translate("MainWindow",  "Rc_Kalman"))
        self.SelectChannel_combo.setItemText(5, _translate("MainWindow",  "Zc_Kalman"))
        self.SelectChannel_combo.setItemText(6, _translate("MainWindow",  "Rc_Mirnov"))
        self.SelectChannel_combo.setItemText(7, _translate("MainWindow",  "Zc_Mirnov"))
        self.SelectChannel_combo.setItemText(8, _translate("MainWindow",  "Rc_Probes"))
        self.SelectChannel_combo.setItemText(9, _translate("MainWindow",  "Zc_Probes"))
               
        self.SelectChannel_combo_2.setItemText(0, _translate("MainWindow", "Vertical Current"))
        self.SelectChannel_combo_2.setItemText(1, _translate("MainWindow", "Plasma Current"))
        self.SelectChannel_combo_2.setItemText(2, _translate("MainWindow", "Primary Current"))
        self.SelectChannel_combo_2.setItemText(3, _translate("MainWindow", "Horizontal Current"))
        self.SelectChannel_combo_2.setItemText(4, _translate("MainWindow",  "Rc_Kalman"))
        self.SelectChannel_combo_2.setItemText(5, _translate("MainWindow",  "Zc_Kalman"))
        self.SelectChannel_combo_2.setItemText(6, _translate("MainWindow",  "Rc_Mirnov"))
        self.SelectChannel_combo_2.setItemText(7, _translate("MainWindow",  "Zc_Mirnov"))
        self.SelectChannel_combo_2.setItemText(8, _translate("MainWindow",  "Rc_Probes"))
        self.SelectChannel_combo_2.setItemText(9, _translate("MainWindow",  "Zc_Probes"))
        
        self.SelectChannel_combo_3.setItemText(0, _translate("MainWindow", "Horizontal Current"))
        self.SelectChannel_combo_3.setItemText(1, _translate("MainWindow", "Vertical Current"))
        self.SelectChannel_combo_3.setItemText(2, _translate("MainWindow", "Plasma Current"))
        self.SelectChannel_combo_3.setItemText(3, _translate("MainWindow", "Primary Current"))
        self.SelectChannel_combo_3.setItemText(4, _translate("MainWindow",  "Rc_Kalman"))
        self.SelectChannel_combo_3.setItemText(5, _translate("MainWindow",  "Zc_Kalman"))
        self.SelectChannel_combo_3.setItemText(6, _translate("MainWindow",  "Rc_Mirnov"))
        self.SelectChannel_combo_3.setItemText(7, _translate("MainWindow",  "Zc_Mirnov"))
        self.SelectChannel_combo_3.setItemText(8, _translate("MainWindow",  "Rc_Probes"))
        self.SelectChannel_combo_3.setItemText(9, _translate("MainWindow",  "Zc_Probes"))
        
        self.SelectChannel_combo_4.setItemText(0, _translate("MainWindow", "Primary Current"))
        self.SelectChannel_combo_4.setItemText(1, _translate("MainWindow", "Vertical Current"))
        self.SelectChannel_combo_4.setItemText(2, _translate("MainWindow", "Horizontal Current"))
        self.SelectChannel_combo_4.setItemText(3, _translate("MainWindow", "Plasma Current"))
        self.SelectChannel_combo_4.setItemText(4, _translate("MainWindow",  "Rc_Kalman"))
        self.SelectChannel_combo_4.setItemText(5, _translate("MainWindow",  "Zc_Kalman"))
        self.SelectChannel_combo_4.setItemText(6, _translate("MainWindow",  "Rc_Mirnov"))
        self.SelectChannel_combo_4.setItemText(7, _translate("MainWindow",  "Zc_Mirnov"))
        self.SelectChannel_combo_4.setItemText(8, _translate("MainWindow",  "Rc_Probes"))
        self.SelectChannel_combo_4.setItemText(9, _translate("MainWindow",  "Zc_Probes"))

        self.comboBox_6.setItemText(0, _translate("MainWindow",  "Rc_Kalman"))
        self.comboBox_6.setItemText(1, _translate("MainWindow",  "Zc_Kalman"))
        self.comboBox_6.setItemText(2, _translate("MainWindow",  "Rc_Mirnov"))
        self.comboBox_6.setItemText(3, _translate("MainWindow",  "Zc_Mirnov"))
        self.comboBox_6.setItemText(4, _translate("MainWindow",  "Rc_Probes"))
        self.comboBox_6.setItemText(5, _translate("MainWindow",  "Zc_Probes"))
        self.comboBox_6.setItemText(6, _translate("MainWindow", "Primary Current"))
        self.comboBox_6.setItemText(7, _translate("MainWindow", "Vertical Current"))
        self.comboBox_6.setItemText(8, _translate("MainWindow", "Horizontal Current"))
        self.comboBox_6.setItemText(9, _translate("MainWindow", "Plasma Current"))
 

        self.comboBox.setItemText(0, _translate("MainWindow",  "Rc_Mirnov"))
        self.comboBox.setItemText(1, _translate("MainWindow",  "Rc_Kalman"))
        self.comboBox.setItemText(2, _translate("MainWindow",  "Zc_Kalman"))
        self.comboBox.setItemText(3, _translate("MainWindow",  "Zc_Mirnov"))
        self.comboBox.setItemText(4, _translate("MainWindow",  "Rc_Probes"))
        self.comboBox.setItemText(5, _translate("MainWindow",  "Zc_Probes"))
        self.comboBox.setItemText(6, _translate("MainWindow", "Primary Current"))
        self.comboBox.setItemText(7, _translate("MainWindow", "Vertical Current"))
        self.comboBox.setItemText(8, _translate("MainWindow", "Horizontal Current"))
        self.comboBox.setItemText(9, _translate("MainWindow", "Plasma Current"))
        
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Database), _translate("MainWindow", "Database"))
        self.Start_Button.setText(_translate("MainWindow", "START"))
        self.Restart_button.setText(_translate("MainWindow", "RESTART"))
        self.ShotNumber_combo.setItemText(0, _translate("MainWindow", "MAX_EVENT"))
        self.ShotNumber_combo.setItemText(1, _translate("MainWindow", "45995"))
        self.ShotNumber_combo.setItemText(2, _translate("MainWindow", "45994"))
        self.ShotNumber_combo.setItemText(3, _translate("MainWindow", "45993"))
        self.ShotNumber_combo.setItemText(4, _translate("MainWindow", "45992"))
        self.ShotNumber_combo.setItemText(5, _translate("MainWindow", "46133"))
        self.ShotNumber_combo.setItemText(6, _translate("MainWindow", "45971"))
        self.ShotNumber_combo.setItemText(7, _translate("MainWindow", "45860"))
        self.ShotNumber_combo.setItemText(8, _translate("MainWindow", "-1"))
        self.ShotNumber_combo.setItemText(9, _translate("MainWindow", "-50"))
        self.ShotNumber_combo.setItemText(10, _translate("MainWindow", "-100"))

    def start_on_button_clicked(self):
        import numpy as np
        #check max event
        if(self.ShotNumber_combo.currentText() == 'MAX_EVENT'):
            SHOT_NUMBER = 0
        else:
            SHOT_NUMBER = self.ShotNumber_combo.currentText()
        try:
            shot_number = int(SHOT_NUMBER)
            data_mirnv_corr, data_mirnv_corr_flux, Ip_magn_corr_value, times, sumIfil_database, R_fromProbes_value, z_fromProbes_value, R_fromMirnov_value, z_fromMirnov_value, MAX_EXENT, tomo_value = getDataFromDatabase(int(shot_number))
            if(self.ShotNumber_combo.currentText() == 'MAX_EVENT'):
                self.textBrowser.append("\nGetting data from the MAX_EVENT")
                self.textBrowser.append("\nMAX_EVENT = " + str(MAX_EXENT))
            elif(int(self.ShotNumber_combo.currentText()) < 0):
                self.textBrowser.append("\nGetting data from the (MAX_EVENT"+self.ShotNumber_combo.currentText()+")")
                self.textBrowser.append("\nSHOT = " + str(MAX_EXENT))
            else:
                self.textBrowser.append("\nGetting data from the shot number: " + self.ShotNumber_combo.currentText())      
            self.textBrowser.append("\nGetting data from database... ")
            B_Mirnov = np.divide((data_mirnv_corr_flux),(49*50*1e-6))
            np.save("B_mirnov.npy", B_Mirnov)
            np.save("I_plasma.npy", Ip_magn_corr_value)
            times_sec = 1e-6 * times
            plotterRM = self.R_mirnov
            item_RM = self.R_mirnov.plot(times_sec, R_fromMirnov_value, pen='b')
            np.save("time_out.npy", times)
            np.save("R_fromMirnov.npy",R_fromMirnov_value)
             
            #mirnov
            plotterRM.setXRange(np.min(times_sec), np.max(times_sec))
            plotterRM.setYRange(np.min(R_fromMirnov_value), np.max(R_fromMirnov_value))
            plotterRM.setLabel('left', 'R', units='m')
            plotterRM.setLabel('bottom', 'Time', units='s')
            plotterRM.showGrid(x=True,y=True,alpha = 2.0)
            plotterRM.addItem(item_RM)
            
            plotterzM = self.Z_mirnov
            item_zM = self.Z_mirnov.plot(times_sec, z_fromMirnov_value, pen='r')
            plotterzM.setXRange(np.min(times_sec), np.max(times_sec))
            plotterzM.setYRange(np.min(z_fromMirnov_value), np.max(z_fromMirnov_value))
            plotterzM.setLabel('left', 'Z', units='m')
            plotterzM.setLabel('bottom', 'Time', units='s')
            plotterzM.showGrid(x=True,y=True,alpha = 2.0)
            np.save("z_fromMirnov.npy",z_fromMirnov_value)
            plotterzM.addItem(item_zM)
            #€lectric probes
            plotterREP = self.R_probes
            item_REP = self.R_probes.plot(times_sec, R_fromProbes_value, pen='b')
            np.save("R_fromProbes.npy", R_fromProbes_value)
            plotterREP.setXRange(np.min(times_sec), np.max(times_sec))
            plotterREP.setYRange(np.min(R_fromProbes_value), np.max(R_fromProbes_value))
            plotterREP.setLabel('left', 'R', units='m')
            plotterREP.setLabel('bottom', 'Time', units='s')
            plotterREP.showGrid(x=True,y=True,alpha = 2.0)            
            plotterREP.addItem(item_REP)
            
            plotterzEP = self.Z_probes
            item_zEP = self.Z_probes.plot(times_sec, z_fromProbes_value, pen='r')
            np.save("z_fromProbes.npy", z_fromProbes_value)
            plotterzEP.setXRange(np.min(times_sec), np.max(times_sec))
            plotterzEP.setYRange(np.min(z_fromProbes_value), np.max(z_fromProbes_value))
            plotterzEP.setLabel('left', 'Z', units='m')
            plotterzEP.setLabel('bottom', 'Time', units='s')
            plotterzEP.showGrid(x=True,y=True,alpha = 2.0) 
            plotterzEP.addItem(item_zEP)
            #tomography
            #R_fromTomgraphy_value, z_fromTomgraphy_value = tomo_centroid(tomo_value)
            R_fromTomgraphy_value = np.load("R_fromProbes.npy")
            z_fromTomgraphy_value = np.load("z_fromProbes.npy")
            plotterRTom = self.R_tomgraphy
            item_RTom = self.R_tomgraphy.plot(times_sec, R_fromTomgraphy_value, pen='b')
            np.save("R_fromTomography.npy", R_fromTomgraphy_value)
            plotterRTom.setXRange(np.min(times_sec), np.max(times_sec))
            plotterRTom.setYRange(np.min(R_fromTomgraphy_value), np.max(R_fromTomgraphy_value))
            plotterRTom.setLabel('left', 'R', units='m')
            plotterRTom.setLabel('bottom', 'Time', units='s')
            plotterRTom.showGrid(x=True,y=True,alpha = 2.0) 
            plotterRTom.addItem(item_RTom)
            
            plotterzTom = self.Z_tomography
            item_zTom = self.Z_tomography.plot(times_sec, z_fromTomgraphy_value, pen='r')
            np.save("z_fromTomography.npy", z_fromTomgraphy_value)
            plotterzTom.setXRange(np.min(times_sec), np.max(times_sec))
            plotterzTom.setYRange(np.min(z_fromTomgraphy_value), np.max(z_fromTomgraphy_value))
            plotterzTom.setLabel('left', 'Z', units='m')
            plotterzTom.setLabel('bottom', 'Time', units='s')
            plotterzTom.showGrid(x=True,y=True,alpha = 2.0)
            plotterzTom.addItem(item_zTom)            
        #ANIMATION OF CENTROID
            self.RZ_mirnov = DynamicMirnovCanvas(self.Centroid_position, width=5, height=4, dpi=55) # QtWidgets.QGraphicsView(self.Centroid_position)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(self.RZ_mirnov.sizePolicy().hasHeightForWidth())
            self.RZ_mirnov.setSizePolicy(sizePolicy)
            self.RZ_mirnov.setObjectName("RZ_mirnov")
            self.verticalLayout_3.addWidget(self.RZ_mirnov)
            self.RZ_probes = DynamicProbesCanvas(self.Centroid_position, width=5, height=4, dpi=55) # QtWidgets.QGraphicsView(self.Centroid_position)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(self.RZ_probes.sizePolicy().hasHeightForWidth())
            self.RZ_probes.setSizePolicy(sizePolicy)
            self.RZ_probes.setObjectName("RZ_probes")
            self.verticalLayout_3.addWidget(self.RZ_probes)
            self.RZ_tomography = DynamicTomographyCanvas(self.Centroid_position, width=5, height=4, dpi=55) # QtWidgets.QGraphicsView(self.Centroid_position)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(self.RZ_tomography.sizePolicy().hasHeightForWidth())
            self.RZ_tomography.setSizePolicy(sizePolicy)
            self.RZ_tomography.setObjectName("RZ_tomography")
            self.verticalLayout_3.addWidget(self.RZ_tomography)
            #END ANIMATION
            #MIRNOV VS CENTROID POSITION
            self.graphicsView_5 = DynamicMinovVsPostionCanvas(self.tab)#QtWidgets.QGraphicsView(self.tab)
            self.graphicsView_5.setObjectName("graphicsView_5")
            self.gridLayout_4.addWidget(self.graphicsView_5, 0, 0, 1, 1)
            #END
            self.Restart_button.clicked.connect(lambda: self.restart_on_button_clicked())
        except ValueError:
            print("You must to write a good shot number!")
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error!")
            msg.setInformativeText("You must to write a good hot number!(exp. 45993, 46093)")
            msg.setWindowTitle("Invalid ShotNumber")
            msg.exec_()

    def restart_on_button_clicked(self):
        self.textBrowser.append("Restarting...")
        python = sys.executable
        os.execl(python, python, * sys.argv)

    def plot_on_button_clicked(self):
        #check max event
        if(self.ShotNumber_combo.currentText() == 'MAX_EVENT'):
            SHOT_NUMBER = 0
        else:
            SHOT_NUMBER = self.ShotNumber_combo.currentText()
        try:
            ShotNumber = int(SHOT_NUMBER)
            Name_1 = str(self.SelectChannel_combo.currentText())
            Name_2 = str(self.SelectChannel_combo_2.currentText())
            Name_3 = str(self.SelectChannel_combo_3.currentText())
            Name_4 = str(self.SelectChannel_combo_4.currentText())
            Channel_1 = select_channel(Name_1)
            Channel_2 = select_channel(Name_2)
            Channel_3 = select_channel(Name_3)
            Channel_4 = select_channel(Name_4)
            self.textBrowser.append("\nPlot1, selected channel: " + Channel_1)
            self.textBrowser.append("\nPlot2, selected channel: " + Channel_2)
            self.textBrowser.append("\nPlot3, selected channel: " + Channel_3)
            self.textBrowser.append("\nPlot4, selected channel: " + Channel_4)
            
            if Channel_1 == "NameError":
                print("error name plot 1")
            else:
                print(Channel_1)
                data_1, dt_1 = getDataFromChannel(Channel_1, ShotNumber)
                time_1 = 1e-6 * dt_1
                plotter_1 = self.graphicsView
                plotter_1.setLabel('bottom', 'Time', units='s')
                item_1 = self.graphicsView.plot(time_1, data_1, pen='r', width = 3)
                plotter_1.setXRange(np.min(time_1), np.max(time_1))
                plotter_1.setYRange(np.min(data_1), np.max(data_1))
                plotter_1.showGrid(x=True,y=True,alpha = 2.0)
                plotter_1.addItem(item_1)
    
            if Channel_2 == "NameError":
                print("error name plot 2")
            else:
                print(Channel_2)
                data_2, dt_2 = getDataFromChannel(Channel_2, ShotNumber)
                time_2 = 1e-6 * dt_2
                plotter_2 = self.graphicsView_2
                plotter_2.setLabel('bottom', 'Time', units='s')
                item_2 = self.graphicsView_2.plot(time_2, data_2, pen='b', width = 3)
                plotter_2.setXRange(np.min(time_2), np.max(time_2))
                plotter_2.setYRange(np.min(data_2), np.max(data_2))
                plotter_2.showGrid(x=True,y=True,alpha = 2.0)
                plotter_2.addItem(item_2)
    
            if Channel_3 == "NameError":
                print("error name plot 3")
            else:
                print(Channel_3)
                data_3, dt_3 = getDataFromChannel(Channel_3, ShotNumber)
                time_3 = 1e-6 * dt_3
                plotter_3 = self.graphicsView_3
                plotter_3.setLabel('bottom', 'Time', units='s')
                item_3 = self.graphicsView_3.plot(time_3, data_3, pen='g', width = 3)
                plotter_3.setXRange(np.min(time_3), np.max(time_3))
                plotter_3.setYRange(np.min(data_3), np.max(data_3))
                plotter_3.showGrid(x=True,y=True,alpha = 2.0)
                plotter_3.addItem(item_3)
    
            if Channel_4 == "NameError":
                print("error name plot 4")
            else:
                print(Channel_4)
                data_4, dt_4 = getDataFromChannel(Channel_4, ShotNumber)
                time_4 = 1e-6 * dt_4
                plotter_4 = self.graphicsView_4
                plotter_4.setLabel('bottom', 'Time', units='s')
                item_4 = self.graphicsView_4.plot(time_4, data_4, pen='m', width = 3)
                plotter_4.setXRange(np.min(time_4), np.max(time_4))
                plotter_4.setYRange(np.min(data_4), np.max(data_4))
                plotter_4.showGrid(x=True,y=True,alpha = 2.0)
                plotter_4.addItem(item_4)

            self.Clear_button.clicked.connect(lambda: self.Clear_on_button_clicked(plotter_1, item_1, plotter_2, item_2, plotter_3, item_3, plotter_4, item_4))
        except ValueError:
            print("You must to write a good shot number!")
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error!")
            msg.setInformativeText("You must to write a good hot number!(exp. 45993, 46093)")
            msg.setWindowTitle("Invalid ShotNumber")
            msg.exec_()
            
    def Clear_on_button_clicked(self, plotter_1, item_1, plotter_2, item_2, plotter_3, item_3, plotter_4, item_4):
        self.textBrowser.append("\nCleaning the plot...")
        plotter_1.removeItem(item_1)
        plotter_1.setXRange(0, 1e+06)
        plotter_1.setYRange(-1, 1)

        plotter_2.removeItem(item_2)
        plotter_2.setXRange(0, 1e+06)
        plotter_2.setYRange(-1, 1)

        plotter_3.removeItem(item_3)
        plotter_3.setXRange(0, 1e+06)
        plotter_3.setYRange(-1, 1)

        plotter_4.removeItem(item_4)
        plotter_4.setXRange(0, 1e+06)
        plotter_4.setYRange(-1, 1)

    def plot_compare_on_button_clicked(self):
        #check max event
        if(self.ShotNumber_combo.currentText() == 'MAX_EVENT'):
            SHOT_NUMBER = 0
        else:
            SHOT_NUMBER = self.ShotNumber_combo.currentText()
        try:
            ShotNumber = int(SHOT_NUMBER)
            Name_1 = str(self.comboBox.currentText())
            Name_2 = str(self.comboBox_6.currentText())
            Channel_1 = select_channel(Name_1)
            Channel_2 = select_channel(Name_2)
            self.textBrowser.append("\nPlot1, selected channel: " + Channel_1)
            self.textBrowser.append("\nPlot2, selected channel: " + Channel_2)
            
            if Channel_1 == "NameError" and Channel_2 == "NameError":
                print("error name plot 1 and error name plot 2")
            else:
                print(Channel_1)
                data_1, dt_1 = getDataFromChannel(Channel_1, ShotNumber)
                time_1 = 1e-6 * dt_1
                plotter_1 = self.graphicsView_6
                plotter_1.setLabel('bottom', 'Time', units='s')
                item_1 = self.graphicsView_6.plot(time_1, data_1, pen='r', width = 3)
                plotter_1.setXRange(np.min(time_1), np.max(time_1))
                plotter_1.setYRange(np.min(data_1), np.max(data_1))
                plotter_1.showGrid(x=True,y=True,alpha = 2.0)
                plotter_1.addItem(item_1)
                
                print(Channel_2)
                data_2, dt_2 = getDataFromChannel(Channel_2, ShotNumber)
                time_2 = 1e-6 * dt_2
                plotter_2 = self.graphicsView_7
                plotter_2.setLabel('bottom', 'Time', units='s')
                item_2 = self.graphicsView_7.plot(time_2, data_2, pen='b', width = 3)
                plotter_2.setXRange(np.min(time_2), np.max(time_2))
                plotter_2.setYRange(np.min(data_2), np.max(data_2))
                plotter_2.showGrid(x=True,y=True,alpha = 2.0)
                plotter_2.addItem(item_2)
                
                item_3 = self.graphicsView_8.plot(time_1, data_1, pen='r', width = 3, name=Name_1)
                item_4 = self.graphicsView_8.plot(time_2, data_2, pen='b', width = 3, name=Name_2)
                plotter_3 = self.graphicsView_8
                plotter_3.setLabel('bottom', 'Time', units='s')
                plotter_3.setXRange(np.min(time_2), np.max(time_2))
                plotter_3.setYRange(np.min(data_2), np.max(data_2))
                plotter_3.showGrid(x=True,y=True,alpha = 2.0)
                #plotter_3.addLegend()
                plotter_3.addItem(item_3, Name_1)
                plotter_3.addItem(item_4, Name_2)


            self.Clear_button_Compare.clicked.connect(lambda: self.Clear_Compare_on_button_clicked(plotter_1, item_1, plotter_2, item_2, plotter_3, item_3, item_4))
        except ValueError:
            print("You must to write a good shot number!")
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error!")
            msg.setInformativeText("You must to write a good hot number!(exp. 45993, 46093)")
            msg.setWindowTitle("Invalid ShotNumber")
            msg.exec_()
    def Clear_Compare_on_button_clicked(self, plotter_1, item_1, plotter_2, item_2, plotter_3, item_3, item_4):
        self.textBrowser.append("\nCleaning the compare plot...")
        plotter_1.removeItem(item_1)
        plotter_1.setXRange(0, 1e+06)
        plotter_1.setYRange(-1, 1)

        plotter_2.removeItem(item_2)
        plotter_2.setXRange(0, 1e+06)
        plotter_2.setYRange(-1, 1)

        plotter_3.removeItem(item_3)
        plotter_3.removeItem(item_4)
        plotter_3.clear()
        plotter_3.setXRange(0, 1e+06)
        plotter_3.setYRange(-1, 1)

           
def main():
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    ui.Start_Button.clicked.connect(lambda: ui.start_on_button_clicked())
    ui.Plot_Button.clicked.connect(lambda: ui.plot_on_button_clicked())
    ui.Plot_Button_Compare.clicked.connect(lambda: ui.plot_compare_on_button_clicked())
    MainWindow.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()