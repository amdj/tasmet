#!/usr/bin/python
import sys
from os.path import isfile
sys.path.append('.')
from PyQt4 import QtCore, QtGui
from tasmet_main_ui import Ui_TASMET
from .. import *
from tasmet_modelclass import TaSMETModel

class StartQT4(QtGui.QMainWindow):
    def __init__(self, parent=None):
        setTASMETTracer(15)
        QtGui.QWidget.__init__(self, parent)
        self.filename=''
        self.model=TaSMETModel()
        self.saved=False
        
        self.ui = Ui_TASMET()
        self.ui.setupUi(self)
        self.ui.gas.addItems([ gas.title() for gas in self.model.gases])

        self.connectRed(self.ui.Nf)
        self.connectRed(self.ui.freq) 
        self.connectRed(self.ui.T0)
        self.connectRed(self.ui.p0) 
        self.connectRed(self.ui.m0)
        # self.connect(self,QtCore.SIGNAL('triggered()'),self.exit)
        self.setGcGreen()
        self.ui.updateGc.clicked.connect(self.updateModel)
        self.ui.file_open.triggered.connect(self.loadModel)
        self.ui.file_save.triggered.connect(self.saveModel)
        self.ui.file_exit.triggered.connect(self.closeEvent)
        self.ui.file_saveas.triggered.connect(self.saveAsModel)                        
        self.ui.m0.setEnabled(False)

        # Update user interface
        self.updateUi()
        self.updateWindowTitle()
        
    def connectRed(self,lineEdit):
        lineEdit.textChanged.connect(self.setGcRed)
        # QtCore.QObject.connect(lineEdit,QtCore.SIGNAL("textChanged(Qstring)"),\
                               # self.setGcRed)
    def setGcGreen(self):
        self.ui.gcbox.setStyleSheet("background-color: #CCFFCC;")

    def setGcRed(self):
        self.ui.gcbox.setStyleSheet("background-color: pink;")
        self.saved=False
        self.updateWindowTitle()
        
    def updateWindowTitle(self):
        if self.filename=='':
            filename='New file.tas'
        else:
            filename=self.filename
        windowtitle="TaSMET - " + filename
        if not self.saved:
            windowtitle=windowtitle+' *'
        self.setWindowTitle(windowtitle)
        
    def updateUi(self):
        self.ui.Nf.setText(str(self.model.Nf))
        self.ui.freq.setText(str(self.model.freq))
        self.ui.T0.setText(str(self.model.T0))
        self.ui.p0.setText(str(self.model.p0))
        index = self.ui.gas.findText(self.model.gas, QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.ui.gas.setCurrentIndex(index)        
        self.setGcGreen()
        
    def loadModel(self):
        # Check if file is saved
        if not self.saved:
            msg='File contains unsaved changes. Continue anyway?'
            reply = QtGui.QMessageBox.question(self, 'Message', \
                                               msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:            
                return

        # Open new file
        filename= QtGui.QFileDialog.getOpenFileName(self,'Open existing TaSMET Model','','*.tas')
        if isfile(filename):
           self.filename=filename
           self.model=TaSMETModel.load(filename)
           self.updateUi()
           self.saved=True
           self.updateWindowTitle()
       
    def saveAsModel(self):
        self.updateModel()
        filename=QtGui.QFileDialog.getSaveFileName(self,'Save TaSMET Model','',"TaSMET Model file (*.tas)")
        if filename[-4:].lower()!='.tas':
           filename=filename+'.tas'
        self.model.save(filename)
        self.filename=filename
        self.updateUi()
        self.saved=True
        self.updateWindowTitle()

    def closeEvent(self,event):
        quit_msg = "File %s contains unsaved changes. Save file first?"
        reply = QtGui.QMessageBox.question(self, 'Message', \
                                            quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No,QtGui.QMessageBox.Cancel)
        if reply == QtGui.QMessageBox.Cancel:
            return
        if reply == QtGui.QMessageBox.Yes:
            self.saveModel()
        # QtCore.QCoreApplication.instance().quit()
        self.deleteLater()

    def saveModel(self):
        self.updateModel()
        if not self.saved:
            if self.filename=='':
                return self.saveAsModel()
            self.model.save(self.filename)
            self.saved=True
            self.updateWindowTitle()
            
    def updateModel(self):
        Nf = self.ui.Nf.text()        
        try:
            Nf = int(Nf)
            if Nf<0 or Nf>30:
                raise ValueError()
        except Exception:
            QtGui.QMessageBox.about(self, 'Error','Input Nf can only be an integer between 0 and 30')
            return
        freq = self.ui.freq.text()        
        try:
            freq = float(freq)
        except Exception:
            QtGui.QMessageBox.about(self, 'Error','Input freq can only be a positive real number')
            return
        
        T0 = self.ui.T0.text()        
        try:
            T0 = float(T0)
            if T0<=0 or T0>1200:
                raise ValueError()
        except Exception:
            QtGui.QMessageBox.about(self, 'Error','Input T0 should be 0<T0<=1200 K.')
            return
        p0 = self.ui.p0.text()        
        try:
            p0 = float(p0)
            if p0<=0 or p0>1e10:
                raise ValueError()
        except Exception:
            QtGui.QMessageBox.about(self, 'Error','Input p0 should be 0<p0<=10 GPa.')
            return
        # gas=self.ui.
        # self.gc=Globalconf(Nf,freq,
        self.setGcGreen()
        self.model.Nf=Nf
        self.model.freq=freq
        self.model.gas=self.ui.gas.currentText().lower()
        self.model.T0=T0
        self.model.p0=p0

def runTaSMETGUI():
    app = QtGui.QApplication(sys.argv)
    myapp = StartQT4()
    myapp.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    run()
