from PyQt4 import QtCore, QtGui
from extremalCurve import Ui_Extremal_Curve
from delayed_filter import DelayedFilter
from base_dialog_widget import BaseDialogWidget
from calpha_choose_chain_to_load_form import CAlphaChooseChainToLoadForm
from seq_model.Chain import Chain
from shutil import copyfile
from os.path import splitext

#import scipy.optimize as optimize
import numpy as np

#import EMAN2.libpyEM.EMAN2
#from EMAN2.libpyEM.EMAN2 import EMAN2db
#from EMAN2.programs import e2version
import os

class ExtremeCurveForm(BaseDialogWidget):
    def __init__(self, main, volumeViewer, parent=None):
        BaseDialogWidget.__init__(self, main, 
                                  "&Extremal Curve", 
                                  "Apply extremal curve skeletonization to density map", 
                                  "perform_VolumeBinarySkeletonization", 
                                  "actions-volume-skeletonization-extremal", 
                                  "actions-volume-skeletonization", 
                                  False, parent)
        self.app = main
        self.viewer = volumeViewer
        self.maxCurveHashes = []
        self.minCurveHashes = []
        self.saddleCurveHashes = []
        self.maxCurve1 = []
        self.maxCurve2 = []
        self.minCurve1 = []
        self.minCurve2 = []
        self.saddleCurve1 = []
        self.saddleCurve2 = []
        self.oldHideCurve = []
        self.skeleton = self.viewer.renderer.getVolume()
        self.vertexPositions = []
        self.quadIndices = []
        self.maxPointSaliencies = []
        self.maxPointIntensities = []
        self.maxPointEigenvalues = []
        self.minPointSaliencies = []
        self.minPointIntensities = []
        self.minPointEigenvalues = []
        self.saddlePointSaliencies = []
        self.saddlePointIntensities = []
        self.saddlePointEigenvalues = []

        self.pointSaliencyMax = -1000.0
        self.pointSaliencyMin = 1000.0
        self.curveSaliencyMax = -1000.0
        self.curveSaliencyMin = 1000.0
        self.surfaceSaliencyMax = -1000.0
        self.surfaceSaliencyMin = 1000.0
        self.minIntensity = 1000.0
        self.maxIntensity = -1000.0
        self.minEigenvalue = 1000.0
        self.maxEigenvalue = -1000.0
        self.fileName = ""
        self.minI = 0;
        self.maxI = 0;
        self.minE = 0;
        self.maxE = 0;

        #self.connect(self.viewer, QtCore.SIGNAL("modelLoaded()"), self.modelLoaded)
        #self.connect(self.viewer, QtCore.SIGNAL("modelUnloaded()"), self.modelUnloaded)
        self.createUI()
        
        self.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.findExtremalCurve)
        self.connect(self.ui.checkBox_22, QtCore.SIGNAL("clicked()"), self.maxCurve)
        self.connect(self.ui.checkBox_24, QtCore.SIGNAL("clicked()"), self.minCurve)
        self.connect(self.ui.checkBox_23, QtCore.SIGNAL("clicked()"), self.saddleCurve)
        #self.connect(self.ui.checkBox_19, QtCore.SIGNAL("clicked()"), self.hideCurve)
        self.connect(self.ui.horizontalSlider, QtCore.SIGNAL("valueChanged(int)"), self.hideCurveAndSurfaceUpdate)
        self.connect(self.ui.horizontalSlider_2, QtCore.SIGNAL("valueChanged(int)"), self.hideCurveAndSurfaceUpdate)
        self.connect(self.ui.horizontalSlider_3, QtCore.SIGNAL("valueChanged(int)"), self.hideCurveAndSurfaceUpdate)
        self.connect(self.ui.horizontalSlider_4, QtCore.SIGNAL("valueChanged(int)"), self.hideCurveAndSurfaceUpdate)
        self.connect(self.ui.horizontalSlider_5, QtCore.SIGNAL("valueChanged(int)"), self.hideCurveAndSurfaceUpdate)
        self.connect(self.ui.horizontalSlider_6, QtCore.SIGNAL("valueChanged(int)"), self.hideCurveAndSurfaceUpdate)
        #self.connect(self.ui.horizontalSlider_7, QtCore.SIGNAL("valueChanged(int)"), self.updateCoilHeight)
        self.connect(self.ui.horizontalSlider_8, QtCore.SIGNAL("valueChanged(int)"), self.updateCoilHeight)

        self.connect(self.ui.checkBox, QtCore.SIGNAL("clicked()"), self.hideCurveAndSurfaceUpdate)
        self.connect(self.ui.checkBox_2, QtCore.SIGNAL("clicked()"), self.hideCurveAndSurfaceUpdate)
        self.connect(self.ui.checkBox_3, QtCore.SIGNAL("clicked()"), self.hideCurveAndSurfaceUpdate)

        #self.connect(self.ui.doubleSpinBox_5, QtCore.SIGNAL("valueChanged(double)"), self.hideSurface)
        #self.connect(self.ui.doubleSpinBox_6, QtCore.SIGNAL("valueChanged(double)"), self.hideSurface)
        #self.connect(self.ui.doubleSpinBox_7, QtCore.SIGNAL("valueChanged(double)"), self.hideSurface)
        #self.connect(self.ui.doubleSpinBox_8, QtCore.SIGNAL("valueChanged(double)"), self.hideSurface)
        #self.connect(self.ui.doubleSpinBox_9, QtCore.SIGNAL("valueChanged(double)"), self.hideSurface)
        #self.connect(self.ui.doubleSpinBox_10, QtCore.SIGNAL("valueChanged(double)"), self.hideSurface)

        self.connect(self.ui.checkBox_25, QtCore.SIGNAL("clicked()"), self.displayMinCheck)
        self.connect(self.ui.checkBox_26, QtCore.SIGNAL("clicked()"), self.displayMaxCheck)
        #self.connect(self.ui.checkBox_20, QtCore.SIGNAL("clicked()"), self.hideSurface)

        self.connect(self.ui.checkBox_21, QtCore.SIGNAL("clicked()"), self.getMinPoints)
        self.connect(self.ui.checkBox_27, QtCore.SIGNAL("clicked()"), self.getMaxPoints)
        self.connect(self.ui.checkBox_4, QtCore.SIGNAL("clicked()"), self.getSaddlePoints)
        self.connect(self.ui.pushButton_2, QtCore.SIGNAL("clicked()"), self.openSkeleton)
        #self.connect(self.ui.pushButton_3, QtCore.SIGNAL("clicked()"), self.writeMaxCurve)
        self.connect(self.ui.pushButton_5, QtCore.SIGNAL("clicked()"), self.writeMaxCurve)
        #self.connect(self.ui.pushButton_4, QtCore.SIGNAL("clicked()"), self.maxCurveHelices)
        #self.connect(self.ui.checkBox_5, QtCore.SIGNAL("clicked()"), self.setHelixDisplay)
        #self.connect(self.ui.checkBox_6, QtCore.SIGNAL("clicked()"), self.setHelixDisplayPoints)
        #self.connect(self.ui.doubleSpinBox_2, QtCore.SIGNAL("valueChanged(double)"), self.setHelixThreshold)

        self.connect(self.ui.checkBox_8, QtCore.SIGNAL("clicked()"), self.setHelixDisplay)
        self.connect(self.ui.checkBox_9, QtCore.SIGNAL("clicked()"), self.setHelixDisplayPoints)
        self.connect(self.ui.doubleSpinBox_4, QtCore.SIGNAL("valueChanged(double)"), self.setHelixThreshold)

        self.connect(self.ui.checkBox_7, QtCore.SIGNAL("clicked()"), self.setEllipsoidVisibility)

        self.connect(self.ui.doubleSpinBox_3, QtCore.SIGNAL("valueChanged(double)"), self.setEllipsoidScale)
        #self.connect(self.ui.spinBox_2, QtCore.SIGNAL("valueChanged(int)"), self.changeSegmentThreshold )
        self.connect(self.ui.spinBox_3, QtCore.SIGNAL("valueChanged(int)"), self.changeSegmentThreshold )
        self.connect(self.ui.horizontalSlider_8, QtCore.SIGNAL("valueChanged(int)"), self.changeOrthoMin )
        self.connect(self.ui.horizontalSlider_9, QtCore.SIGNAL("valueChanged(int)"), self.changeOrthoMax )
        self.connect(self.ui.doubleSpinBox_5, QtCore.SIGNAL("valueChanged(double)"), self.updateCoilHeight )

        #self.connect(self.ui.checkBox_22, QtCore.SIGNAL("clicked()"), self.maxCurve)
        #self.createActions()

    def createUI(self):
      self.ui = Ui_Extremal_Curve()
      self.ui.setupUi(self)

    def changeOrthoMin(self):
      outOf1 = float(self.ui.horizontalSlider_8.value()) / 100.0
      self.app.viewers['calpha'].renderer.setOrthoDistMin(outOf1)
      self.app.viewers['calpha'].emitModelChanged()

    def changeOrthoMax(self):
      outOf1 = float(self.ui.horizontalSlider_9.value()) / 100.0
      self.app.viewers['calpha'].renderer.setOrthoDistMax(outOf1)
      self.app.viewers['calpha'].emitModelChanged()


    def setHelixDisplay(self):
      if self.ui.checkBox_8.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setHelixDisplay(True)
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setHelixDisplay(False)
        self.app.viewers['calpha'].emitModelChanged()

    def updateCoilHeight(self):
      self.app.viewers['calpha'].renderer.setHelixCoilHeight(self.ui.doubleSpinBox_5.value())
      self.app.viewers['calpha'].emitModelChanged()

    def setEllipsoidVisibility(self):
      if self.ui.checkBox_7.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setEllipsoidDisplay(True)
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setEllipsoidDisplay(False)
        self.app.viewers['calpha'].emitModelChanged()

    def setHelixDisplayPoints(self):
      if self.ui.checkBox_9.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setHelixDisplayPoints(True)
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setHelixDisplayPoints(False)
        self.app.viewers['calpha'].emitModelChanged()

    def setHelixThreshold(self):
      self.app.viewers['calpha'].renderer.setHelixDisplayThreshold(float(self.ui.doubleSpinBox_4.value()))
      self.app.viewers['calpha'].emitModelChanged()

    def setEllipsoidScale(self):
      self.app.viewers['calpha'].renderer.setScaleEllipsoid(float(self.ui.doubleSpinBox_3.value()) )
      self.app.viewers['calpha'].emitModelChanged()

    def generateAtoms(self, fileName):
      def setupChain(mychain):            
            self.app.viewers['calpha'].main_chain = mychain
            self.app.viewers['calpha'].loadedChains.append(mychain)
            mychain.setViewer(self.app.viewers['calpha'])
            #mychain.addCalphaBondsPathwalker()
            #mychain.addSideChainBonds()
            renderer = self.app.viewers['calpha'].renderer
            for i in mychain.unsortedResidueRange():
                for atomName in mychain[i].getAtomNames():
                    atom = mychain[i].getAtom(atomName)
                    if atom:
                        atom = renderer.addAtom(atom)
                        mychain[i].addAtomObject(atom)
      
      self.atomFileName = fileName
      fileNameTemp = self.atomFileName
      self.app.viewers['calpha'].whichChainID = None
      filename = unicode(self.atomFileName)
      if filename.split('.')[-1].lower() == 'pdb':
          #dlg = CAlphaChooseChainToLoadForm(unicode(self.atomFileName))
          if True:
              self.app.viewers['calpha'].whichChainID = 'ALL'
              #self.app.viewers['calpha'].whichChainID = dlg.whichChainID
              #if not self.atomFileName.isEmpty():
              if(self.app.viewers['calpha'].loaded):
                  self.app.viewers['calpha'].unloadData()
              
              self.atomFileName = fileNameTemp
                    
              if self.app.viewers['calpha'].whichChainID == 'ALL':
                  mychainKeys = Chain.loadAllChainsPathwalker(str(self.atomFileName), qparent=self.app)
                  #mychainKeys = Chain.loadAllChains(str(self.atomFileName), qparent=self.app)
                  for chainKey in mychainKeys:
                    setupChain(Chain.getChain(chainKey))
              else:
                  mychain = Chain.__loadFromPDBPathwalker(str(self.atomFileName), qparent=self.app, whichChainID = self.app.viewers['calpha'].whichChainID)
                  #mychain = Chain.load(str(self.atomFileName), qparent=self.app, whichChainID = self.app.viewers['calpha'].whichChainID)
                  setupChain(mychain)
        
              if not self.app.viewers['calpha'].loaded:
                  #self.app.viewers['calpha'].setDisplayStyle(6)
                  #self.app.viewers['calpha'].displayStyle = 6
                  #self.app.viewers['calpha'].renderer.setDisplayStyle(6)
                  #self.app.viewers['calpha'].setAtomColorsAndVisibility(6)
                  #self.app.viewers['calpha'].modelChangedPathwalker()
                  self.app.viewers['calpha'].dirty = False
                  self.app.viewers['calpha'].loaded = True
                  self.app.viewers['calpha'].setAtomColorsAndVisibility(self.app.viewers['calpha'].displayStyle)                        
                  self.app.viewers['calpha'].emitModelLoadedPreDraw()
                  #self.app.viewers['calpha'].emitModelPathwalker()
                  self.app.viewers['calpha'].emitModelLoaded()
                  self.app.viewers['calpha'].emitViewerSetCenter()

    def displayMinCheck(self):
      if self.ui.checkBox_25.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setMinSurfaceOn(True)
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setMinSurfaceOn(False)
        self.app.viewers['calpha'].emitModelChanged()

    def displayMaxCheck(self):
      if self.ui.checkBox_26.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setMaxSurfaceOn(True)
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setMaxSurfaceOn(False)
        self.app.viewers['calpha'].emitModelChanged()

    def hideSurface(self):
      #if self.ui.checkBox_20.checkState() == QtCore.Qt.Checked:
      #  self.app.viewers['calpha'].renderer.setExtremalSurfaceHide(True)
      #else:
      #  self.app.viewers['calpha'].renderer.setExtremalSurfaceHide(False)
      #self.app.viewers['calpha'].renderer.setExtremalParams(float(self.ui.horizontalSlider.value()), float(self.ui.horizontalSlider_2.value()), float(self.ui.horizontalSlider_3.value()), float(self.ui.horizontalSlider_4.value()), float(self.ui.horizontalSlider_6.value()), float(self.ui.horizontalSlider_5.value()))
      self.updateParams()
      saliencyChecked = False
      if self.ui.checkBox.checkState() == QtCore.Qt.Checked:
        saliencyChecked = True
      intensityChecked = False
      if self.ui.checkBox_2.checkState() == QtCore.Qt.Checked:
        intensityChecked = True
      eigenChecked = False
      if self.ui.checkBox_3.checkState() == QtCore.Qt.Checked:
        eigenChecked = True
      self.app.viewers['calpha'].renderer.setExtremalChecks(saliencyChecked, intensityChecked, eigenChecked)
      self.app.viewers['calpha'].emitModelChanged()     

    def findExtremalCurve(self):
      
      self.skeleton = self.viewer.renderer.performExtremalCurve2016(int(self.ui.spinBox.value()), 1, 1, 1, 1, float(self.ui.doubleSpinBox.value()), 1.0, 3);
      
      #self.skeleton = skeleton
      self.app.viewers["skeleton"].loadVolume(self.skeleton)
      self.app.viewers["skeleton"].visualizationOptions.ui.checkBoxModelVisible.setChecked(False)
      self.app.viewers["skeleton"].visualizationOptions.ui.checkBoxModel2Visible.setChecked(False)
      self.app.viewers["skeleton"].visualizationOptions.ui.checkBoxModel3Visible.setChecked(False)
      self.app.viewers['calpha'].visualizationOptions.ui.checkBoxShowAtoms.setChecked(False)
      #self.ui.checkBox_25.setChecked(True)
      #self.ui.checkBox_26.setChecked(True)
      #self.displaySurface()

      self.generateAtoms("extremal.pdb")

      self.maxCurveHashes = self.skeleton.getMaxHashes()
      self.minCurveHashes = self.skeleton.getMinHashes()
      self.saddleCurveHashes = self.skeleton.getSaddleHashes()
      self.updateParams()
      
      self.vertexPositions = self.skeleton.getVertexPos()
      self.sendRendererPositions()
      
      self.maxCurve1 = self.skeleton.getExtremalBonds1()
      self.maxCurve2 = self.skeleton.getExtremalBonds2()

      self.maxCurveSaliencies = self.skeleton.getMaxSaliencies()
      self.maxCurveIntensities = self.skeleton.getMaxIntensities()
      self.maxCurveEigenvalues = self.skeleton.getMaxEigenvalues()
      #for i in self.maxCurveEigenvalues:
      #  print "mce " + str(i)
      self.maxCurveEigenvalues0 = self.skeleton.getMaxCurveEigenvalues0()
      #for i in self.maxCurveEigenvalues0:
      #  print "mce0 " + str(i)
      self.maxCurveEigenvalues1 = self.skeleton.getMaxCurveEigenvalues1()
      #for i in self.maxCurveEigenvalues1:
      #  print "mce1 " + str(i)
      #maxMaxCurveEigenvalue2 = max(self.maxCurveEigenvalues)
      #maxMaxCurveEigenvalue0 = max(self.maxCurveEigenvalues0)
      #maxMaxCurveEigenvalue1 = max(self.maxCurveEigenvalues1)
      #minMaxCurveEigenvalue2 = min(self.maxCurveEigenvalues)
      #minMaxCurveEigenvalue0 = min(self.maxCurveEigenvalues0)
      #minMaxCurveEigenvalue1 = min(self.maxCurveEigenvalues1)
      maxMaxCurveEigenvalue2 = -1000.0
      maxMaxCurveEigenvalue0 = -1000.0
      maxMaxCurveEigenvalue1 = -1000.0
      minMaxCurveEigenvalue2 = 1000.0
      minMaxCurveEigenvalue0 = 1000.0
      minMaxCurveEigenvalue1 = 1000.0
      for i in range(len(self.maxCurveEigenvalues)):
        if self.maxCurveEigenvalues0[i] > maxMaxCurveEigenvalue0 and self.maxCurveEigenvalues0[i] > 0:
          maxMaxCurveEigenvalue0 = self.maxCurveEigenvalues0[i]
        if self.maxCurveEigenvalues1[i] > maxMaxCurveEigenvalue1 and self.maxCurveEigenvalues1[i] > 0:
          maxMaxCurveEigenvalue1 = self.maxCurveEigenvalues1[i]
        if self.maxCurveEigenvalues[i] > maxMaxCurveEigenvalue2 and self.maxCurveEigenvalues[i] > 0:
          maxMaxCurveEigenvalue2 = self.maxCurveEigenvalues[i]

        if self.maxCurveEigenvalues0[i] < minMaxCurveEigenvalue0 and self.maxCurveEigenvalues0[i] > 0:
          minMaxCurveEigenvalue0 = self.maxCurveEigenvalues0[i]
        if self.maxCurveEigenvalues1[i] < minMaxCurveEigenvalue1 and self.maxCurveEigenvalues1[i] > 0:
          minMaxCurveEigenvalue1 = self.maxCurveEigenvalues1[i]
        if self.maxCurveEigenvalues[i] < minMaxCurveEigenvalue2 and self.maxCurveEigenvalues[i] > 0:
          minMaxCurveEigenvalue2 = self.maxCurveEigenvalues[i]
      #print str(maxMaxCurveEigenvalue2) + " " + str(maxMaxCurveEigenvalue0) + " " + str(maxMaxCurveEigenvalue1) + " " + str(minMaxCurveEigenvalue2) + " " + str(minMaxCurveEigenvalue0) + " " + str(minMaxCurveEigenvalue1)

      self.app.viewers['calpha'].renderer.setMaxCurveEigenMaxesAndMins(maxMaxCurveEigenvalue0, maxMaxCurveEigenvalue1, maxMaxCurveEigenvalue2, minMaxCurveEigenvalue0, minMaxCurveEigenvalue1, minMaxCurveEigenvalue2)

      self.maxCurveEigenvectors = self.skeleton.getMaxCurveEigenvectors()
      self.minCurveEigenvectors = self.skeleton.getMinCurveEigenvectors()
      self.saddleCurveEigenvectors = self.skeleton.getSaddleCurveEigenvectors()


      self.minCurve1 = self.skeleton.getMinCurveBonds1()
      self.minCurve2 = self.skeleton.getMinCurveBonds2()

      self.minCurveSaliencies = self.skeleton.getMinSaliencies()
      self.minCurveIntensities = self.skeleton.getMinIntensities()
      self.minCurveEigenvalues = self.skeleton.getMinEigenvalues()
      self.minCurveEigenvalues0 = self.skeleton.getMinCurveEigenvalues0()
      self.minCurveEigenvalues1 = self.skeleton.getMinCurveEigenvalues1()

      self.saddleCurve1 = self.skeleton.getSaddleCurveBonds1()
      self.saddleCurve2 = self.skeleton.getSaddleCurveBonds2()

      self.saddleCurveSaliencies = self.skeleton.getSaddleSaliencies()
      self.saddleCurveIntensities = self.skeleton.getSaddleIntensities()
      self.saddleCurveEigenvalues = self.skeleton.getSaddleEigenvalues()
      self.saddleCurveEigenvalues0 = self.skeleton.getSaddleCurveEigenvalues0()
      self.saddleCurveEigenvalues1 = self.skeleton.getSaddleCurveEigenvalues1()
      self.app.viewers['calpha'].renderer.setScaleEllipsoid(float(self.ui.doubleSpinBox_3.value()) )

      # Test to store the curve out
      #outfile=open('testCurve.off','w')

      self.app.viewers['calpha'].renderer.setExtremalMode(True)
      print "storing max curve..."
      self.findMaxCurve()
      print "storing min curve..."
      self.findMinCurve()
      print "storing saddle curve..."
      self.findSaddleCurve()
      self.findMinMaxEigenvalues()
      #self.app.viewers['calpha'].renderer.setExtremalParams(float(self.ui.doubleSpinBox_5.value()), float(self.ui.doubleSpinBox_6.value()), float(self.ui.doubleSpinBox_7.value()), float(self.ui.doubleSpinBox_8.value()), float(self.ui.doubleSpinBox_9.value()), float(self.ui.doubleSpinBox_10.value()))      
      
      #self.app.viewers['calpha'].renderer.setExtremalParams(float(self.ui.horizontalSlider.value()), float(self.ui.horizontalSlider_2.value()), float(self.ui.horizontalSlider_3.value()), float(self.ui.horizontalSlider_4.value()), float(self.ui.horizontalSlider_6.value()), float(self.ui.horizontalSlider_5.value()))
      self.app.viewers['calpha'].renderer.sortMaxCurveDimensions() 

      print "storing surfaces..."
      self.quadSaliencies = self.skeleton.getQuadSaliencies()
      self.quadIntensities = self.skeleton.getQuadIntensities()
      self.quadEigenvalues = self.skeleton.getQuadEigenvalues()
      self.quadIndices = self.skeleton.getQuadIndices()
      self.displayQuads = self.skeleton.getDisplayQuads()
      self.displayNormals = self.skeleton.getExtremalNormals()
      self.quadTypes = self.skeleton.getQuadTypes()
      self.quadEigenvalues0 = self.skeleton.getQuadEigenvalues0()
      self.quadEigenvalues1 = self.skeleton.getQuadEigenvalues1()
      self.quadEigenvectors = self.skeleton.getQuadEigenvectors()
      self.findQuads()
      print "storing min points..."
      self.extremalMinPoints = self.skeleton.getExtremalMinPoints()
      self.minPointSaliencies = self.skeleton.getMinPointSaliencies()
      self.minPointIntensities = self.skeleton.getMinPointIntensities()
      self.minPointEigenvalues = self.skeleton.getMinPointEigenvalues()
      self.minPointEigenvalues0 = self.skeleton.getMinPointEigenvalues0()
      self.minPointEigenvalues1 = self.skeleton.getMinPointEigenvalues1()
      self.minPointEigenvectors = self.skeleton.getMinPointEigenvectors()
      for i in range(len(self.extremalMinPoints)):
          self.app.viewers['calpha'].renderer.makeMin(self.extremalMinPoints[i], True, self.minPointSaliencies[3*i], self.minPointSaliencies[3*i+1], self.minPointSaliencies[3*i+2], self.minPointIntensities[i], self.minPointEigenvalues[i], self.minPointEigenvalues0[i], self.minPointEigenvalues1[i])
          #self.app.viewers['calpha'].renderer.addMinPtEigenVector(self.extremalMinPoints[i], self.minPointEigenvectors[9*i], self.minPointEigenvectors[9*i+1], self.minPointEigenvectors[9*i+2], self.minPointEigenvectors[9*i+3], self.minPointEigenvectors[9*i+4], self.minPointEigenvectors[9*i+5], self.minPointEigenvectors[9*i+6], self.minPointEigenvectors[9*i+7], self.minPointEigenvectors[9*i+8])

      print "storing max points..."
      self.extremalMaxPoints = self.skeleton.getExtremalMaxPoints()
      self.maxPointSaliencies = self.skeleton.getMaxPointSaliencies()
      self.maxPointIntensities = self.skeleton.getMaxPointIntensities()
      self.maxPointEigenvalues = self.skeleton.getMaxPointEigenvalues()
      self.maxPointEigenvalues0 = self.skeleton.getMaxPointEigenvalues0()
      self.maxPointEigenvalues1 = self.skeleton.getMaxPointEigenvalues1()
      self.maxPointEigenvectors = self.skeleton.getMaxPointEigenvectors()

      print str(len(self.extremalMaxPoints)) + " " + str(len(self.maxPointSaliencies)) + " " + str(len(self.maxPointIntensities)) + " " + str(len(self.maxPointEigenvalues))
      for i in range(len(self.extremalMaxPoints)):
          self.app.viewers['calpha'].renderer.makeMax(self.extremalMaxPoints[i], True, self.maxPointSaliencies[3*i], self.maxPointSaliencies[3*i+1], self.maxPointSaliencies[3*i+2], self.maxPointIntensities[i], self.maxPointEigenvalues[i], self.maxPointEigenvalues0[i], self.maxPointEigenvalues1[i])
          #self.app.viewers['calpha'].renderer.addMaxPtEigenVector(self.extremalMaxPoints[i], self.maxPointEigenvectors[9*i], self.maxPointEigenvectors[9*i+1], self.maxPointEigenvectors[9*i+2], self.maxPointEigenvectors[9*i+3], self.maxPointEigenvectors[9*i+4], self.maxPointEigenvectors[9*i+5], self.maxPointEigenvectors[9*i+6], self.maxPointEigenvectors[9*i+7], self.maxPointEigenvectors[9*i+8])
      print "storing saddle points..."
      self.extremalSaddlePoints = self.skeleton.getExtremalSaddlePoints()
      self.saddlePointSaliencies = self.skeleton.getSaddlePointSaliencies()
      self.saddlePointIntensities = self.skeleton.getSaddlePointIntensities()
      self.saddlePointEigenvalues = self.skeleton.getSaddlePointEigenvalues()
      self.saddlePointEigenvalues0 = self.skeleton.getSaddlePointEigenvalues0()
      self.saddlePointEigenvalues1 = self.skeleton.getSaddlePointEigenvalues1()
      self.saddlePointEigenvectors = self.skeleton.getSaddlePointEigenvectors()
      for i in range(len(self.extremalSaddlePoints)):
          #print str(len(self.extremalSaddlePoints)) + " " + str(len(self.saddlePointSaliencies)) + " " + str(len(self.saddlePointIntensities)) + " " + str(len(self.saddlePointEigenvalues))
          self.app.viewers['calpha'].renderer.makeSaddle(self.extremalSaddlePoints[i], True, self.saddlePointSaliencies[3*i], self.saddlePointSaliencies[3*i+1], self.saddlePointSaliencies[3*i+2], self.saddlePointIntensities[i], self.saddlePointEigenvalues[i], self.saddlePointEigenvalues0[i], self.saddlePointEigenvalues1[i])
          #self.app.viewers['calpha'].renderer.addSaddlePtEigenVector(self.extremalSaddlePoints[i], self.saddlePointEigenvectors[9*i], self.saddlePointEigenvectors[9*i+1], self.saddlePointEigenvectors[9*i+2], self.saddlePointEigenvectors[9*i+3], self.saddlePointEigenvectors[9*i+4], self.saddlePointEigenvectors[9*i+5], self.saddlePointEigenvectors[9*i+6], self.saddlePointEigenvectors[9*i+7], self.saddlePointEigenvectors[9*i+8])

      self.app.viewers['calpha'].renderer.findHelices()
      #self.setMinMaxHideParameters()
      self.updateParams()
      self.app.viewers['calpha'].renderer.setHelixDisplayThreshold(float(self.ui.doubleSpinBox_4.value()))
      self.app.viewers['calpha'].renderer.setHelixCoilHeight(self.ui.doubleSpinBox_5.value())
      self.app.viewers['calpha'].renderer.setSegmentThreshold(int(self.ui.spinBox_3.value()))
      self.app.viewers['calpha'].renderer.setOrthoDistMin(float(self.ui.horizontalSlider_8.value()) / 100.0)
      self.app.viewers['calpha'].renderer.setOrthoDistMax(float(self.ui.horizontalSlider_9.value()) / 100.0)
      self.app.viewers['calpha'].renderer.averageEllipsoidScale()
      print "Save skeleton.txt ..."
      skeletonFileName = QtGui.QFileDialog.getSaveFileName(self, self.tr("Save Data"), "", self.tr('Atom Positions (*.ske)\nAll files (*.*)'))
      skeletonFileName = str(skeletonFileName)
      if skeletonFileName != '':
        copyfile('extremal.pdb', splitext(skeletonFileName)[0]+'.pdb')
        self.saveSkeleton(skeletonFileName)
      print "rendering ..."
      #os.system("python prehelix.py --pdbin extremal.pdb --mapin map.mrc --output path1.pdb")

    def fitFunc(self, data, a, b, c):
      return data[:,0]*a + data[:,0]*b + c

    def writeMaxCurve(self):
      helicesSize = self.app.viewers['calpha'].renderer.findHelices()
      fitCurves = []
      for i in range(helicesSize):
        currentStrand = self.app.viewers['calpha'].renderer.findMaxCurveStrand(i)
        strandTuples = []
        for j in range(len(currentStrand)/3):
          strandTuples.append((currentStrand[3*j], currentStrand[3*j+1], currentStrand[3*j+2]))
        fitCurves.append(strandTuples)
      curveCoeffA = []
      curveCoeffB = []
      curveCoeffC = []
      #curveCoeffD = []
      ords = []
      for strand in fitCurves:
        x = []
        y = []
        x1 = []
        x2 = []
        for dataItem in strand:
          x1.append(dataItem[0])
          x2.append(dataItem[1])
          y.append(dataItem[2])


        x1 = np.array(x1, np.float64);
        x2 = np.array(x2, np.float64);
        X = np.array([x1, x2], np.float64);
        n = np.max(X.shape)    
        X = np.vstack([np.ones(n), X]).T
        y = np.array(y, np.float64)
        data = np.concatenate((x1[:, np.newaxis], x2[:, np.newaxis], y[:, np.newaxis]), axis=1)
        #data = np.array([x1, x2, y])
        datamean = data.mean(axis=0)
        uu, dd, vv = np.linalg.svd(data - datamean)
        print "svd"
        print uu
        print dd
        print vv
        print vv[0]
        print "end svd"
        for coeff in vv[0]:
          ords.append(coeff)



        #X = np.column_stack(x+[[1]*len(x[0])])
        #print X
        #print X.shape
        #print x1.shape
        #print x2.shape
        #print y.shape
        beta_hat = np.linalg.lstsq(X,y)[0]
        
        #print beta_hat

        #A = np.array(strand)
        #guess = (1,1,1)
        #params, pcov = optimize.curve_fit(self.fitFunc, A[:,:2], A[:,2], guess)
        #print params

        curveCoeffA.append(beta_hat[0])
        curveCoeffB.append(beta_hat[1])
        curveCoeffC.append(beta_hat[2])
        #curveCoeffD.append(c)
      self.app.viewers['calpha'].renderer.findCurveMinsMaxes()
      for i in range(helicesSize):
        self.app.viewers['calpha'].renderer.addCurveCoeff(curveCoeffA[i], curveCoeffB[i], curveCoeffC[i])
        self.app.viewers['calpha'].renderer.addOrd(ords[3*i], ords[3*i+1], ords[3*i+2])
        #self.app.viewers['calpha'].renderer.addCurveCoeff(curveCoeffA[i], curveCoeffB[i], curveCoeffC[i], curveCoeffD[i])
      self.app.viewers['calpha'].renderer.interpolateHelixPoints()
      #self.app.viewers['calpha'].renderer.normalizedHelicesMinsAndMaxes()
      self.app.viewers['calpha'].emitModelChanged()

    def changeSegmentThreshold(self):
      self.app.viewers['calpha'].renderer.setSegmentThreshold(int(self.ui.spinBox_3.value()))
      self.app.viewers['calpha'].emitModelChanged()

    def maxCurveHelices(self):
      #numHelices = self.app.viewers['calpha'].renderer.findMaxHelices()
      os.system("python prehelix.py --pdbin maxCurves.pdb --mapin map.mrc --output skeletonH.pdb")
      #for i in range(numHelices):
      #os.system("python prehelix.py --pdbin tempFiles/" + str(8) + ".pdb --mapin map.mrc --output tempFiles/helix" + str(8) + ".pdb")

    def findQuads(self):
      for i in range(len(self.quadIndices)/4):
        #print "quad eigenvalues " + str(self.quadEigenvalues0[i]) + " " + str(self.quadEigenvalues1[i]) + " " + str(self.quadEigenvalues[i])
        self.app.viewers['calpha'].renderer.addQuadSurface(self.quadIndices[4*i], self.quadIndices[4*i+1], self.quadIndices[4*i+2], self.quadIndices[4*i+3], self.displayNormals[3*i], self.displayNormals[3*i+1], self.displayNormals[3*i+2], self.quadSaliencies[3*i], self.quadSaliencies[3*i+1], self.quadSaliencies[3*i+2], self.quadIntensities[i], self.quadEigenvalues[i], self.quadTypes[i], self.quadEigenvalues0[i], self.quadEigenvalues1[i])
        self.app.viewers['calpha'].renderer.addQuadEigenvectors(self.quadIndices[4*i], self.quadIndices[4*i+1], self.quadIndices[4*i+2], self.quadIndices[4*i+3], self.quadEigenvectors[9*i], self.quadEigenvectors[9*i+1], self.quadEigenvectors[9*i+2], self.quadEigenvectors[9*i+3], self.quadEigenvectors[9*i+4], self.quadEigenvectors[9*i+5], self.quadEigenvectors[9*i+6], self.quadEigenvectors[9*i+7], self.quadEigenvectors[9*i+8])
      self.app.viewers['calpha'].emitModelChanged()
      #for i in range(len(self.quadIndices)/4):
       # self.app.viewers['calpha'].renderer.addQuadEigenvectors(self.quadIndices[4*i], self.quadIndices[4*i+1], self.quadIndices[4*i+2], self.quadIndices[4*i+3], self.quadEigenvectors[9*i], self.quadEigenvectors[9*i+1], self.quadEigenvectors[9*i+2], self.quadEigenvectors[9*i+3], self.quadEigenvectors[9*i+4], self.quadEigenvectors[9*i+5], self.quadEigenvectors[9*i+6], self.quadEigenvectors[9*i+7], self.quadEigenvectors[9*i+8])
        #self.app.viewers['calpha'].renderer.addQuadSurface(self.quadIndices[4*i], self.quadIndices[4*i+1], self.quadIndices[4*i+2],self.quadIndices[4*i+3], self.displayQuads[4*i], self.displayQuads[4*i+1], self.displayQuads[4*i+2], self.displayQuads[4*i+3], self.displayNormals[3*i], self.displayNormals[3*i+1], self.displayNormals[3*i+2], self.quadSaliencies[3*i], self.quadSaliencies[3*i+1], self.quadSaliencies[3*i+2], self.quadIntensities[i], self.quadEigenvalues[i], self.quadTypes[i])

        #self.app.viewers['calpha'].renderer.addQuadSurface(self.quadIndices[4*i], self.quadIndices[4*i+1], self.quadIndices[4*i+2],self.quadIndices[4*i+3], self.displayQuads[4*i], self.displayQuads[4*i+1], self.displayQuads[4*i+2], self.displayQuads[4*i+3], self.displayNormals[3*i], self.displayNormals[3*i+1], self.displayNormals[3*i+2], self.quadSaliencies[3*i], self.quadSaliencies[3*i+1], self.quadSaliencies[3*i+2], self.quadIntensities[i], self.quadEigenvalues[i], self.quadTypes[i])
      #if self.ui.checkBox_22.checkState() == QtCore.Qt.Checked:
       # for i in range(len(bonds1)):
        #  print "bonds " + str(bonds1[i]) + " " + str(bonds2[i])
                    #print str(bonds1[i]) + " " + str(bonds2[i])
         # self.app.viewers['calpha'].renderer.drawAddedBond(bonds1[i], bonds2[i])
          #self.app.viewers['calpha'].emitModelChanged()

    def setMinMaxHideParameters(self):
      maxCurvePointSaliencies = []
      maxCurveCurveSaliencies = []
      maxCurveSurfaceSaliencies = []
      for i in range(len(self.maxCurveSaliencies)/3):
        maxCurvePointSaliencies.append(self.maxCurveSaliencies[3*i])
        maxCurveCurveSaliencies.append(self.maxCurveSaliencies[3*i+1])
        maxCurveSurfaceSaliencies.append(self.maxCurveSaliencies[3*i+2])
      minCurvePointSaliencies = []
      minCurveCurveSaliencies = []
      minCurveSurfaceSaliencies = []
      for i in range(len(self.minCurveSaliencies)/3):
        minCurvePointSaliencies.append(self.minCurveSaliencies[3*i+2])
        minCurveCurveSaliencies.append(self.minCurveSaliencies[3*i+1])
        minCurveSurfaceSaliencies.append(self.minCurveSaliencies[3*i])
      saddleCurvePointSaliencies = []
      saddleCurveCurveSaliencies = []
      saddleCurveSurfaceSaliencies = []
      for i in range(len(self.saddleCurveSaliencies)/3):
        saddleCurvePointSaliencies.append(self.saddleCurveSaliencies[3*i+2])
        saddleCurveCurveSaliencies.append(self.saddleCurveSaliencies[3*i+1])
        saddleCurveSurfaceSaliencies.append(self.saddleCurveSaliencies[3*i])
      minPointPointSaliencies = []
      minPointCurveSaliencies = []
      minPointSurfaceSaliencies = []
      for i in range(len(self.minPointSaliencies)/3):
        minPointPointSaliencies.append(self.minPointSaliencies[3*i+2])
        minPointCurveSaliencies.append(self.minPointSaliencies[3*i+1])
        minPointSurfaceSaliencies.append(self.minPointSaliencies[3*i])
      maxPointPointSaliencies = []
      maxPointCurveSaliencies = []
      maxPointSurfaceSaliencies = []
      for i in range(len(self.maxPointSaliencies)/3):
        maxPointPointSaliencies.append(self.maxPointSaliencies[3*i+2])
        maxPointCurveSaliencies.append(self.maxPointSaliencies[3*i+1])
        maxPointSurfaceSaliencies.append(self.maxPointSaliencies[3*i])
      saddlePointPointSaliencies = []
      saddlePointCurveSaliencies = []
      saddlePointSurfaceSaliencies = []
      for i in range(len(self.saddlePointSaliencies)/3):
        saddlePointPointSaliencies.append(self.saddlePointSaliencies[3*i+2])
        saddlePointCurveSaliencies.append(self.saddlePointSaliencies[3*i+1])
        saddlePointSurfaceSaliencies.append(self.saddlePointSaliencies[3*i])

      quadPointSaliencies = []
      quadCurveSaliencies = []
      quadSurfaceSaliencies = []
      for i in range(len(self.quadSaliencies)/3):
        quadPointSaliencies.append(self.quadSaliencies[3*i+2])
        quadCurveSaliencies.append(self.quadSaliencies[3*i+1])
        quadSurfaceSaliencies.append(self.quadSaliencies[3*i])
      self.pointSaliencyMax = max([max(maxCurvePointSaliencies), max(minCurvePointSaliencies), max(saddleCurvePointSaliencies), max(quadPointSaliencies), max(minPointPointSaliencies), max(maxPointPointSaliencies), max(saddlePointPointSaliencies)])
      self.pointSaliencyMin = min([min(maxCurvePointSaliencies), min(minCurvePointSaliencies), min(saddleCurvePointSaliencies), min(quadPointSaliencies), min(minPointPointSaliencies), min(maxPointPointSaliencies), min(saddlePointPointSaliencies)])
      self.curveSaliencyMax = max([max(maxCurveCurveSaliencies), max(minCurveCurveSaliencies), max(saddleCurveCurveSaliencies), max(quadCurveSaliencies), max(minPointCurveSaliencies), max(maxPointCurveSaliencies), max(saddlePointCurveSaliencies)])
      self.curveSaliencyMin = min([min(maxCurveCurveSaliencies), min(minCurveCurveSaliencies), min(saddleCurveCurveSaliencies), min(quadCurveSaliencies), min(minPointCurveSaliencies), min(maxPointCurveSaliencies), min(saddlePointCurveSaliencies)])
      self.surfaceSaliencyMax = max([max(maxCurveSurfaceSaliencies), max(minCurveSurfaceSaliencies), max(saddleCurveSurfaceSaliencies), max(quadSurfaceSaliencies), max(minPointSurfaceSaliencies), max(maxPointSurfaceSaliencies), max(saddlePointSurfaceSaliencies)])
      self.surfaceSaliencyMin = min([min(maxCurveSurfaceSaliencies), min(minCurveSurfaceSaliencies), min(saddleCurveSurfaceSaliencies), min(quadSurfaceSaliencies), min(minPointSurfaceSaliencies), min(maxPointSurfaceSaliencies), min(saddlePointSurfaceSaliencies)])
      self.minIntensity = min(min(self.minCurveIntensities), min(self.maxCurveIntensities), min(self.saddleCurveIntensities), min(self.quadIntensities), min(self.minPointIntensities), min(self.maxPointIntensities), min(self.saddlePointIntensities))
      self.maxIntensity = max(max(self.minCurveIntensities), max(self.maxCurveIntensities), max(self.saddleCurveIntensities), max(self.quadIntensities), max(self.minPointIntensities), max(self.maxPointIntensities), max(self.saddlePointIntensities))
      self.minEigenvalue = min(min(self.minCurveEigenvalues), min(self.maxCurveEigenvalues), min(self.saddleCurveEigenvalues), min(self.quadEigenvalues), min(self.minPointEigenvalues), min(self.maxPointEigenvalues), min(self.saddlePointEigenvalues))
      self.maxEigenvalue = max(max(self.minCurveEigenvalues), max(self.maxCurveEigenvalues), max(self.saddleCurveEigenvalues), max(self.quadEigenvalues), max(self.minPointEigenvalues), max(self.maxPointEigenvalues), max(self.saddlePointEigenvalues))

      #quadPointSaliencyMax = max(quadPointSaliencies)
      #quadCurveSaliencyMax = max(quadCurveSaliencies)
      #quadSurfaceSalienciesMax = max(quadSurfaceSaliencies)
      #quadPointSaliencyMin = min(quadPointSaliencies)
      #quadCurveSaliencyMin = min(quadCurveSaliencies)
      #quadSurfaceSalienciesMin = min(quadSurfaceSaliencies)


    def findSaddleCurve(self):
      print "sizes " + str(len(self.saddleCurveHashes)) + " " + str(len(self.saddleCurve1)) + " " + str(len(self.saddleCurve2)) + " " + str(len(self.saddleCurveSaliencies)) + " " + str(len(self.saddleCurveIntensities)) + " " + str(len(self.saddleCurveEigenvalues))
      for i in range(len(self.saddleCurveHashes)/2):
          self.app.viewers['calpha'].renderer.addSaddleBond(self.saddleCurve1[i], self.saddleCurve2[i], self.saddleCurveHashes[2*i], self.saddleCurveHashes[2*i+1], self.saddleCurveSaliencies[3*i], self.saddleCurveSaliencies[(3*i)+1], self.saddleCurveSaliencies[(3*i)+2], self.saddleCurveIntensities[i], self.saddleCurveEigenvalues[i], self.saddleCurveEigenvalues0[i], self.saddleCurveEigenvalues1[i])
      for i in range(len(self.saddleCurveHashes)/2):
        self.app.viewers['calpha'].renderer.addMinEigenVector(self.saddleCurveEigenvectors[9*i], self.saddleCurveEigenvectors[9*i+1], self.saddleCurveEigenvectors[9*i+2], self.saddleCurveEigenvectors[9*i+3], self.saddleCurveEigenvectors[9*i+4], self.saddleCurveEigenvectors[9*i+5], self.saddleCurveEigenvectors[9*i+6],self.saddleCurveEigenvectors[9*i+7], self.saddleCurveEigenvectors[9*i+8], self.saddleCurveHashes[2*i], self.saddleCurveHashes[2*i+1]);

      self.app.viewers['calpha'].emitModelChanged()

    def sendRendererPositions(self):
      for i in range(len(self.vertexPositions)/3):
        self.app.viewers['calpha'].renderer.addVertexPos(self.vertexPositions[3*i],self.vertexPositions[3*i+1], self.vertexPositions[3*i+2])

    def findMinMaxEigenvalues(self):
      self.minI = self.skeleton.getMinI();
      self.maxI = self.skeleton.getMaxI();
      self.minE = self.skeleton.getMinE();
      self.maxE = self.skeleton.getMaxE();
      self.app.viewers['calpha'].renderer.setMinMaxIntensitiesEigenvalues(self.minI, self.maxI, self.minE, self.maxE)
      self.app.viewers['calpha'].emitModelChanged()

    def saddleCurve(self):
      if self.ui.checkBox_23.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setSaddleOn(True)        
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setSaddleOn(False)
        self.app.viewers['calpha'].emitModelChanged()

    def hideCurveAndSurfaceUpdate(self):
      self.hideCurve()
      self.hideSurface()

    def findMinCurve(self):
      for i in range(len(self.minCurveHashes)/2):
          self.app.viewers['calpha'].renderer.addMinBond(self.minCurve1[i], self.minCurve2[i], self.minCurveHashes[2*i], self.minCurveHashes[2*i+1], self.minCurveSaliencies[3*i], self.minCurveSaliencies[(3*i)+1], self.minCurveSaliencies[(3*i)+2], self.minCurveIntensities[i], self.minCurveEigenvalues[i], self.minCurveEigenvalues0[i], self.minCurveEigenvalues1[i])
      for i in range(len(self.minCurveHashes)/2):
        self.app.viewers['calpha'].renderer.addMinEigenVector(self.minCurveEigenvectors[9*i], self.minCurveEigenvectors[9*i+1], self.minCurveEigenvectors[9*i+2], self.minCurveEigenvectors[9*i+3], self.minCurveEigenvectors[9*i+4], self.minCurveEigenvectors[9*i+5], self.minCurveEigenvectors[9*i+6],self.minCurveEigenvectors[9*i+7], self.minCurveEigenvectors[9*i+8], self.minCurveHashes[2*i], self.minCurveHashes[2*i+1]);

      self.app.viewers['calpha'].emitModelChanged()

    def minCurve(self):
      if self.ui.checkBox_24.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setMinOn(True)        
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setMinOn(False) 
        self.app.viewers['calpha'].emitModelChanged()

    def findMaxCurve(self):
      for i in range(len(self.maxCurveHashes)/2):
        self.app.viewers['calpha'].renderer.addMaxBond(self.maxCurve1[i], self.maxCurve2[i], self.maxCurveHashes[2*i], self.maxCurveHashes[2*i+1], self.maxCurveSaliencies[3*i], self.maxCurveSaliencies[(3*i)+1], self.maxCurveSaliencies[(3*i)+2], self.maxCurveIntensities[i], self.maxCurveEigenvalues0[i], self.maxCurveEigenvalues1[i], self.maxCurveEigenvalues[i])
      for i in range(len(self.maxCurveHashes)/2):
        self.app.viewers['calpha'].renderer.addMaxEigenVector(self.maxCurveEigenvectors[9*i], self.maxCurveEigenvectors[9*i+1], self.maxCurveEigenvectors[9*i+2], self.maxCurveEigenvectors[9*i+3], self.maxCurveEigenvectors[9*i+4], self.maxCurveEigenvectors[9*i+5], self.maxCurveEigenvectors[9*i+6],self.maxCurveEigenvectors[9*i+7], self.maxCurveEigenvectors[9*i+8], self.maxCurveHashes[2*i], self.maxCurveHashes[2*i+1]);
      #for i in range(len(self.maxCurve1)):
          #self.app.viewers['calpha'].renderer.addMaxBond(self.maxCurve1[i], self.maxCurve2[i], self.maxCurveSaliencies[3*i], self.maxCurveSaliencies[(3*i)+1], self.maxCurveSaliencies[(3*i)+2], self.maxCurveIntensities[i], self.maxCurveEigenvalues[i])
      self.app.viewers['calpha'].renderer.findEigenMinMaxes()
      self.app.viewers['calpha'].emitModelChanged()


    def maxCurve(self):
      if self.ui.checkBox_22.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setMaxOn(True)
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setMaxOn(False)
        self.app.viewers['calpha'].emitModelChanged()

    def hideCurve(self):
      #if self.ui.checkBox_19.checkState() == QtCore.Qt.Checked:
       # self.app.viewers['calpha'].renderer.setExtremalHide(True)
      #else:
       # self.app.viewers['calpha'].renderer.setExtremalHide(False)
      #self.app.viewers['calpha'].renderer.setExtremalParams(float(self.ui.horizontalSlider.value()), float(self.ui.horizontalSlider_2.value()), float(self.ui.horizontalSlider_3.value()), float(self.ui.horizontalSlider_4.value()), float(self.ui.horizontalSlider_6.value()), float(self.ui.horizontalSlider_5.value()))
      self.updateParams()
      saliencyChecked = False
      if self.ui.checkBox.checkState() == QtCore.Qt.Checked:
        saliencyChecked = True
      intensityChecked = False
      if self.ui.checkBox_2.checkState() == QtCore.Qt.Checked:
        intensityChecked = True

      eigenChecked = False
      if self.ui.checkBox_3.checkState() == QtCore.Qt.Checked:
        eigenChecked = True
      self.app.viewers['calpha'].renderer.setExtremalChecks(saliencyChecked, intensityChecked, eigenChecked)
      self.app.viewers['calpha'].emitModelChanged()

    def updateParams(self):
      pointRatio = (float(self.ui.horizontalSlider.value())/100.0)
      curveRatio = (float(self.ui.horizontalSlider_2.value())/100.0)
      surfaceRatio = (float(self.ui.horizontalSlider_3.value())/100.0)
      minGeo = (float(self.ui.horizontalSlider_4.value())/100.0)
      maxGeo = (float(self.ui.horizontalSlider_6.value())/100.0)
      eigenValue = (float(self.ui.horizontalSlider_5.value())/100.0)

      #pointRatio = ((float(self.ui.horizontalSlider.value())/100.0)*20.0)-5.0
      #curveRatio = ((float(self.ui.horizontalSlider_2.value())/100.0)*20.0)-5.0
      #surfaceRatio = ((float(self.ui.horizontalSlider_3.value())/100.0)*20.0)-5.0
      #minGeo = ((float(self.ui.horizontalSlider_4.value())/100.0)*20.0)-5.0
      #maxGeo = ((float(self.ui.horizontalSlider_6.value())/100.0)*20.0)-5.0
      #eigenValue = ((float(self.ui.horizontalSlider_5.value())/100.0)*20.0)-5.0
      
      #pointRatio = -1.0 * self.pointSaliencyMin + (3.0*float(self.ui.horizontalSlider.value())/100.0)*(self.pointSaliencyMax - self.pointSaliencyMin)
      #curveRatio = -1.0 * self.curveSaliencyMin + (3.0*float(self.ui.horizontalSlider_2.value())/100.0)*(self.curveSaliencyMax-self.curveSaliencyMin)
      #surfaceRatio = -1.0 * self.surfaceSaliencyMin + (3.0*float(self.ui.horizontalSlider_3.value())/100.0)*(self.surfaceSaliencyMax - self.surfaceSaliencyMin)
      #minGeo = -1.0 * self.minIntensity + (3.0*float(self.ui.horizontalSlider_4.value())/100.0)*(self.maxIntensity - self.minIntensity)
      #maxGeo = -1.0 * self.minIntensity + (3.0*float(self.ui.horizontalSlider_6.value())/100.0)*(self.maxIntensity - self.minIntensity)
      #eigenValue = -1.0 * self.minEigenvalue + (3.0*float(self.ui.horizontalSlider_5.value())/100.0)*(self.maxEigenvalue - self.minEigenvalue)
      self.app.viewers['calpha'].renderer.setExtremalParams(curveRatio, pointRatio, surfaceRatio, minGeo, maxGeo, eigenValue)
      self.app.viewers['calpha'].emitModelChanged()

    def getMinPoints(self):
      if self.ui.checkBox_21.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setMinDisplay(True)
        self.app.viewers['calpha'].emitModelChanged()
      else:
        self.app.viewers['calpha'].renderer.setMinDisplay(False)
        self.app.viewers['calpha'].emitModelChanged()

    def getMaxPoints(self):
      if self.ui.checkBox_27.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setMaxDisplay(True)
        self.app.viewers['calpha'].emitModelChanged()    
      else:
        self.app.viewers['calpha'].renderer.setMaxDisplay(False)
        self.app.viewers['calpha'].emitModelChanged()

    def getSaddlePoints(self):
      if self.ui.checkBox_4.checkState() == QtCore.Qt.Checked:
        self.app.viewers['calpha'].renderer.setSaddleDisplay(True)
        self.app.viewers['calpha'].emitModelChanged()    
      else:
        self.app.viewers['calpha'].renderer.setSaddleDisplay(False)
        self.app.viewers['calpha'].emitModelChanged()

    def openSkeleton(self):
      self.fileName = QtGui.QFileDialog.getOpenFileName(self, self.tr("Open Data"), "", self.tr('ExtCurSke (*.ske)\nAll files (*.*)'))
      if str(self.fileName) == '':
          return

      skeletonData = []
      with open(self.fileName) as f:
        skeletonData = f.readlines()
      self.maxCurveHashes = skeletonData[0].split()
      self.minCurveHashes = skeletonData[1].split()
      self.saddleCurveHashes = skeletonData[2].split()
      self.updateParams()

      pdbFileName = splitext(str(self.fileName))[0] + '.pdb'
      self.generateAtoms(pdbFileName)

      self.vertexPositions = skeletonData[9].split()
      self.vertexPositions = [float(i) for i in skeletonData[9].split()]
      self.sendRendererPositions()

      self.maxCurve1 = skeletonData[3].split()
      self.maxCurve1 = [int(i) for i in self.maxCurve1]

      self.maxCurve2 = skeletonData[4].split()
      self.maxCurve2 = [int(i) for i in self.maxCurve2]

      self.maxCurveSaliencies = skeletonData[31].split()
      self.maxCurveSaliencies = [float(i) for i in self.maxCurveSaliencies]

      self.maxCurveIntensities = skeletonData[32].split()
      self.maxCurveIntensities = [float(i) for i in self.maxCurveIntensities]
      
      self.maxCurveEigenvalues = skeletonData[33].split()
      self.maxCurveEigenvalues = [float(i) for i in self.maxCurveEigenvalues]

      self.minCurve1 = skeletonData[5].split()
      self.minCurve1 = [int(i) for i in self.minCurve1]

      self.minCurve2 = skeletonData[6].split()
      self.minCurve2 = [int(i) for i in self.minCurve2]

      self.minCurveSaliencies = skeletonData[27].split()
      self.minCurveSaliencies = [float(i) for i in self.minCurveSaliencies]

      self.minCurveIntensities = skeletonData[28].split()
      self.minCurveIntensities = [float(i) for i in self.minCurveIntensities]
      
      self.minCurveEigenvalues = skeletonData[29].split()
      self.minCurveEigenvalues = [float(i) for i in self.minCurveEigenvalues]


      self.saddleCurve1 = skeletonData[7].split()
      self.saddleCurve1 = [int(i) for i in self.saddleCurve1]

      self.saddleCurve2 = skeletonData[8].split()
      self.saddleCurve2 = [int(i) for i in self.saddleCurve2]

      self.saddleCurveSaliencies = skeletonData[35].split()
      self.saddleCurveSaliencies = [float(i) for i in self.saddleCurveSaliencies]

      self.saddleCurveIntensities = skeletonData[36].split()
      self.saddleCurveIntensities = [float(i) for i in self.saddleCurveIntensities]

      self.saddleCurveEigenvalues = skeletonData[37].split()
      self.saddleCurveEigenvalues = [float(i) for i in self.saddleCurveEigenvalues]

      self.maxCurveEigenvalues0 = skeletonData[42].split()
      self.maxCurveEigenvalues0 = [float(i) for i in self.maxCurveEigenvalues0]

      self.maxCurveEigenvalues1 = skeletonData[43].split()
      self.maxCurveEigenvalues1 = [float(i) for i in self.maxCurveEigenvalues1]

      self.minCurveEigenvalues0 = skeletonData[44].split()
      self.minCurveEigenvalues0 = [float(i) for i in self.minCurveEigenvalues0]

      self.minCurveEigenvalues1 = skeletonData[45].split()
      self.minCurveEigenvalues1 = [float(i) for i in self.minCurveEigenvalues1]

      self.maxCurveEigenvectors = skeletonData[46].split()
      self.maxCurveEigenvectors = [float(i) for i in self.maxCurveEigenvectors]

      self.minCurveEigenvectors = skeletonData[47].split()
      self.minCurveEigenvectors = [float(i) for i in self.minCurveEigenvectors]

      self.saddlePointEigenvalues0 = skeletonData[48].split()
      self.saddlePointEigenvalues0 = [float(i) for i in self.saddlePointEigenvalues0]

      self.saddlePointEigenvalues1 = skeletonData[49].split()
      self.saddlePointEigenvalues1 = [float(i) for i in self.saddlePointEigenvalues1]

      self.saddleCurveEigenvalues0 = skeletonData[50].split()
      self.saddleCurveEigenvalues0 = [float(i) for i in self.saddleCurveEigenvalues0]

      self.saddleCurveEigenvalues1 = skeletonData[51].split()
      self.saddleCurveEigenvalues1 = [float(i) for i in self.saddleCurveEigenvalues1]

      self.saddleCurveEigenvectors = skeletonData[52].split()
      self.saddleCurveEigenvectors = [float(i) for i in self.saddleCurveEigenvectors]

      self.quadEigenvalues0 = skeletonData[53].split()
      self.quadEigenvalues0 = [float(i) for i in self.quadEigenvalues0]

      self.quadEigenvalues1 = skeletonData[54].split()
      self.quadEigenvalues1 = [float(i) for i in self.quadEigenvalues1]

      self.quadEigenvectors = skeletonData[55].split()
      self.quadEigenvectors = [float(i) for i in self.quadEigenvectors]

      self.minPointEigenvalues0 = skeletonData[56].split()
      self.minPointEigenvalues0 = [float(i) for i in self.minPointEigenvalues0]

      self.minPointEigenvalues1 = skeletonData[57].split()
      self.minPointEigenvalues1 = [float(i) for i in self.minPointEigenvalues1]

      self.maxPointEigenvalues0 = skeletonData[58].split()
      self.maxPointEigenvalues0 = [float(i) for i in self.maxPointEigenvalues0]

      self.maxPointEigenvalues1 = skeletonData[59].split()
      self.maxPointEigenvalues1 = [float(i) for i in self.maxPointEigenvalues1]

      '''
      #self.minPointEigenvectors = skeletonData[60].split()
      #self.minPointEigenvectors = [float(i) for i in self.minPointEigenvectors]

      self.maxPointEigenvectors = skeletonData[61].split()
      self.maxPointEigenvectors = [float(i) for i in self.maxPointEigenvectors]

      self.saddlePointEigenvectors = skeletonData[62].split()
      self.saddlePointEigenvectors = [float(i) for i in self.saddlePointEigenvectors]
      '''

      self.app.viewers['calpha'].renderer.setExtremalMode(True)
      print "finding max curve... "
      self.findMaxCurve()
      print "finding min curve... "
      self.findMinCurve()
      print "finding saddle curve... "
      self.findSaddleCurve()

      self.minI = skeletonData[38].split()
      self.minI = float(self.minI[0])

      self.maxI = skeletonData[39].split()
      self.maxI = float(self.maxI[0])

      self.minE = skeletonData[40].split()
      self.minE = float(self.minE[0])

      self.maxE = skeletonData[41].split()
      self.maxE = float(self.maxE[0])

      self.app.viewers['calpha'].renderer.setMinMaxIntensitiesEigenvalues(self.minI, self.maxI, self.minE, self.maxE)

      #self.findMinMaxEigenvalues()
      print "storing surfaces..."

      self.quadSaliencies = skeletonData[20].split()
      self.quadSaliencies = [float(i) for i in self.quadSaliencies]

      self.quadIntensities = skeletonData[21].split()
      self.quadIntensities = [float(i) for i in self.quadIntensities]

      self.quadEigenvalues = skeletonData[22].split()
      self.quadEigenvalues = [float(i) for i in self.quadEigenvalues]

      self.quadIndices = skeletonData[10].split()
      self.quadIndices = [int(i) for i in self.quadIndices]

      self.displayQuads = skeletonData[23].split()

      self.displayNormals = skeletonData[24].split()
      self.displayNormals = [float(i) for i in self.displayNormals]

      self.quadTypes = skeletonData[25].split()
      self.quadTypes = [int(i) for i in self.quadTypes]

      self.findQuads()

      print "storing min points..."

      self.extremalMinPoints = skeletonData[26].split()
      self.extremalMinPoints = [int(i) for i in self.extremalMinPoints]

      self.minPointSaliencies = skeletonData[14].split()
      self.minPointSaliencies = [float(i) for i in self.minPointSaliencies]

      self.minPointIntensities = skeletonData[15].split()
      self.minPointIntensities = [float(i) for i in self.minPointIntensities]

      self.minPointEigenvalues = skeletonData[16].split()
      self.minPointEigenvalues = [float(i) for i in self.minPointEigenvalues]
      
      for i in range(len(self.extremalMinPoints)):
          self.app.viewers['calpha'].renderer.makeMin(self.extremalMinPoints[i], True, self.minPointSaliencies[3*i], self.minPointSaliencies[3*i+1], self.minPointSaliencies[3*i+2], self.minPointIntensities[i], self.minPointEigenvalues[i], self.minPointEigenvalues0[i], self.minPointEigenvalues1[i])
      
      print "storing max points..."
      self.extremalMaxPoints = skeletonData[30].split()
      self.extremalMaxPoints = [int(i) for i in self.extremalMaxPoints]
      

      self.maxPointSaliencies = skeletonData[11].split()
      self.maxPointSaliencies = [float(i) for i in self.maxPointSaliencies]

      self.maxPointIntensities = skeletonData[12].split()
      self.maxPointIntensities = [float(i) for i in self.maxPointIntensities]

      self.maxPointEigenvalues = skeletonData[13].split()
      self.maxPointEigenvalues = [float(i) for i in self.maxPointEigenvalues]

      for i in range(len(self.extremalMaxPoints)):
          self.app.viewers['calpha'].renderer.makeMax(self.extremalMaxPoints[i], True, self.maxPointSaliencies[3*i], self.maxPointSaliencies[3*i+1], self.maxPointSaliencies[3*i+2], self.maxPointIntensities[i], self.maxPointEigenvalues[i], self.maxPointEigenvalues0[i], self.maxPointEigenvalues1[i])
      
      self.extremalSaddlePoints = skeletonData[34].split()
      self.extremalSaddlePoints = [int(i) for i in self.extremalSaddlePoints]

      self.saddlePointSaliencies = skeletonData[17].split()
      self.saddlePointSaliencies = [float(i) for i in self.saddlePointSaliencies]

      self.saddlePointIntensities = skeletonData[18].split()
      self.saddlePointIntensities = [float(i) for i in self.saddlePointIntensities]

      self.saddlePointEigenvalues = skeletonData[19].split()
      self.saddlePointEigenvalues = [float(i) for i in self.saddlePointEigenvalues]

      print "storing saddle points..."
      for i in range(len(self.extremalSaddlePoints)):
          self.app.viewers['calpha'].renderer.makeSaddle(self.extremalSaddlePoints[i], True, self.saddlePointSaliencies[3*i], self.saddlePointSaliencies[3*i+1], self.saddlePointSaliencies[3*i+2], self.saddlePointIntensities[i], self.saddlePointEigenvalues[i], self.saddlePointEigenvalues0[i], self.saddlePointEigenvalues1[i])

      self.app.viewers['calpha'].renderer.setHelixDisplayThreshold(float(self.ui.doubleSpinBox_4.value()))
      self.app.viewers['calpha'].renderer.setHelixCoilHeight(self.ui.doubleSpinBox_5.value())
      self.app.viewers['calpha'].renderer.setSegmentThreshold(int(self.ui.spinBox_3.value()))
      self.app.viewers['calpha'].renderer.setOrthoDistMin(float(self.ui.horizontalSlider_8.value()) / 100.0)
      self.app.viewers['calpha'].renderer.setOrthoDistMax(float(self.ui.horizontalSlider_9.value()) / 100.0)
      self.app.viewers['calpha'].renderer.setScaleEllipsoid(float(self.ui.doubleSpinBox_3.value()) )
      self.app.viewers['calpha'].renderer.averageEllipsoidScale()

      self.updateParams()
      self.app.viewers['calpha'].emitModelChanged()
      
      

    def saveSkeleton(self, fileName):
      f = open(fileName, 'w')
      for i in range(len(self.maxCurveHashes)):#0
        item = str(self.maxCurveHashes[i])
        if i != len(self.maxCurveHashes)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurveHashes)):#1
        item = str(self.minCurveHashes[i])
        if i != len(self.minCurveHashes)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurveHashes)):#2
        item = str(self.saddleCurveHashes[i])
        if i != len(self.saddleCurveHashes)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxCurve1)):#3
        item = str(self.maxCurve1[i])
        if i != len(self.maxCurve1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxCurve2)):#4
        item = str(self.maxCurve2[i])
        if i != len(self.maxCurve2)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurve1)):#5
        item = str(self.minCurve1[i])
        if i != len(self.minCurve1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurve2)):#6
        item = str(self.minCurve2[i])
        if i != len(self.minCurve2)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurve1)):#7
        item = str(self.saddleCurve1[i])
        if i != len(self.saddleCurve1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurve2)):#8
        item = str(self.saddleCurve2[i])
        if i != len(self.saddleCurve2)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.vertexPositions)):#9
        item = str(self.vertexPositions[i])
        if i != len(self.vertexPositions)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadIndices)):#10
        item = str(self.quadIndices[i])
        if i != len(self.quadIndices)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxPointSaliencies)):#11
        item = str(self.maxPointSaliencies[i])
        if i != len(self.maxPointSaliencies)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)
      
      for i in range(len(self.maxPointIntensities)):#12
        item = str(self.maxPointIntensities[i])
        if i != len(self.maxPointIntensities)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxPointEigenvalues)):#13
        item = str(self.maxPointEigenvalues[i])
        if i != len(self.maxPointEigenvalues)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minPointSaliencies)):#14
        item = str(self.minPointSaliencies[i])
        if i != len(self.minPointSaliencies)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minPointIntensities)):#15
        item = str(self.minPointIntensities[i])
        if i != len(self.minPointIntensities)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minPointEigenvalues)):#16
        item = str(self.minPointEigenvalues[i])
        if i != len(self.minPointEigenvalues)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddlePointSaliencies)):#17
        item = str(self.saddlePointSaliencies[i])
        if i != len(self.saddlePointSaliencies)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddlePointIntensities)):#18
        item = str(self.saddlePointIntensities[i])
        if i != len(self.saddlePointIntensities)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddlePointEigenvalues)):#19
        item = str(self.saddlePointEigenvalues[i])
        if i != len(self.saddlePointEigenvalues)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadSaliencies)):#20
        item = str(self.quadSaliencies[i])
        if i != len(self.quadSaliencies)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadIntensities)):#21
        item = str(self.quadIntensities[i])
        if i != len(self.quadIntensities)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadEigenvalues)):#22
        item = str(self.quadEigenvalues[i])
        if i != len(self.quadEigenvalues)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.displayQuads)):#23
        item = str(self.displayQuads[i])
        if i != len(self.displayQuads)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.displayNormals)):#24
        item = str(self.displayNormals[i])
        if i != len(self.displayNormals)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadTypes)):#25
        item = str(self.quadTypes[i])
        if i != len(self.quadTypes)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.extremalMinPoints)):#26
        item = str(self.extremalMinPoints[i])
        if i != len(self.extremalMinPoints)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurveSaliencies)):#27
        item = str(self.minCurveSaliencies[i])
        if i != len(self.minCurveSaliencies)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurveIntensities)):#28
        item = str(self.minCurveIntensities[i])
        if i != len(self.minCurveIntensities)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurveEigenvalues)):#29
        item = str(self.minCurveEigenvalues[i])
        if i != len(self.minCurveEigenvalues)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.extremalMaxPoints)):#30
        item = str(self.extremalMaxPoints[i])
        if i != len(self.extremalMaxPoints)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxCurveSaliencies)):#31
        item = str(self.maxCurveSaliencies[i])
        if i != len(self.maxCurveSaliencies)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxCurveIntensities)):#32
        item = str(self.maxCurveIntensities[i])
        if i != len(self.maxCurveIntensities)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxCurveEigenvalues)):#33
        item = str(self.maxCurveEigenvalues[i])
        if i != len(self.maxCurveEigenvalues)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.extremalSaddlePoints)):#34
        item = str(self.extremalSaddlePoints[i])
        if i != len(self.extremalSaddlePoints)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurveSaliencies)):#35
        item = str(self.saddleCurveSaliencies[i])
        if i != len(self.saddleCurveSaliencies)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurveIntensities)):#36
        item = str(self.saddleCurveIntensities[i])
        if i != len(self.saddleCurveIntensities)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurveEigenvalues)):#37
        item = str(self.saddleCurveEigenvalues[i])
        if i != len(self.saddleCurveEigenvalues)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      f.write(str(self.minI) + "\n") #38
      f.write(str(self.maxI) + "\n") #39
      f.write(str(self.minE) + "\n") #40
      f.write(str(self.maxE) + "\n") #41

      for i in range(len(self.maxCurveEigenvalues0)):#42
        item = str(self.maxCurveEigenvalues0[i])
        if i != len(self.maxCurveEigenvalues0)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxCurveEigenvalues1)):#43
        item = str(self.maxCurveEigenvalues1[i])
        if i != len(self.maxCurveEigenvalues1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurveEigenvalues0)):#44
        item = str(self.minCurveEigenvalues0[i])
        if i != len(self.minCurveEigenvalues0)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurveEigenvalues1)):#45
        item = str(self.minCurveEigenvalues1[i])
        if i != len(self.minCurveEigenvalues1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxCurveEigenvectors)):#46
        item = str(self.maxCurveEigenvectors[i])
        if i != len(self.maxCurveEigenvectors)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minCurveEigenvectors)):#47
        item = str(self.minCurveEigenvectors[i])
        if i != len(self.minCurveEigenvectors)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddlePointEigenvalues0)):#48
        item = str(self.saddlePointEigenvalues0[i])
        if i != len(self.saddlePointEigenvalues0)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddlePointEigenvalues1)):#49
        item = str(self.saddlePointEigenvalues1[i])
        if i != len(self.saddlePointEigenvalues1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurveEigenvalues0)):#50
        item = str(self.saddleCurveEigenvalues0[i])
        if i != len(self.saddleCurveEigenvalues0)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurveEigenvalues1)):#51
        item = str(self.saddleCurveEigenvalues1[i])
        if i != len(self.saddleCurveEigenvalues1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddleCurveEigenvectors)):#52
        item = str(self.saddleCurveEigenvectors[i])
        if i != len(self.saddleCurveEigenvectors)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadEigenvalues0)):#53
        item = str(self.quadEigenvalues0[i])
        if i != len(self.quadEigenvalues0)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadEigenvalues1)):#54
        item = str(self.quadEigenvalues1[i])
        if i != len(self.quadEigenvalues1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.quadEigenvectors)):#55
        item = str(self.quadEigenvectors[i])
        if i != len(self.quadEigenvectors)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minPointEigenvalues0)):#56
        item = str(self.minPointEigenvalues0[i])
        if i != len(self.minPointEigenvalues0)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.minPointEigenvalues1)):#57
        item = str(self.minPointEigenvalues1[i])
        if i != len(self.minPointEigenvalues1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxPointEigenvalues0)):#58
        item = str(self.maxPointEigenvalues0[i])
        if i != len(self.maxPointEigenvalues0)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxPointEigenvalues1)):#59
        item = str(self.maxPointEigenvalues1[i])
        if i != len(self.maxPointEigenvalues1)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      '''
      for i in range(len(self.minPointEigenvectors)):#60
        item = str(self.minPointEigenvectors[i])
        if i != len(self.minPointEigenvectors)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.maxPointEigenvectors)):#61
        item = str(self.maxPointEigenvectors[i])
        if i != len(self.maxPointEigenvectors)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)

      for i in range(len(self.saddlePointEigenvectors)):#62
        item = str(self.saddlePointEigenvectors[i])
        if i != len(self.saddlePointEigenvectors)-1:
          item += " "
        else:
          item += "\n"
        f.write(item)
      '''

      f.close()


