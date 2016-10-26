# Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
# Author:        Sasakthi S. Abeysinghe (sasakthi@gmail.com)
# Description:   A widget to fit a calpha backbone to the density

# CVS Meta Information: 
#   $Source: /project/mm/cvs/graphics/ssa1/source/Gorgon/src_py/calpha_flexible_fitting_form.py,v $
#   $Revision: 1.17 $
#   $Date: 2012/06/17 13:27:00 $
#   $Author: burrowsd $
#   $State: Exp $
#
# History Log: 
#   $Log: calpha_flexible_fitting_form.py,v $
#   Revision 1.17  2012/06/18 13:27:00  burrowsd
#   Complete re-write of flexible fitting plugin
#
#   Revision 1.16  2010/08/19 23:05:08  chenb
#   Cleaned and commented ribbon diagram code
#
#   Revision 1.15  2010/07/23 18:18:33  heiderp
#   Side chains now transform correctly.  PDB helices now color correctly and rigid initialization bug is fixed
#
#   Revision 1.14  2010/07/22 21:09:07  heiderp
#   Minor updates. Mostly commenting and removing extra material from CurveDeformer.h
#
#   Revision 1.13  2010/07/19 17:29:02  heiderp
#   LARGE update.  Added flexible fitting functionality, lots of logic in FlexibleFittingEngine.h
#
#   Revision 1.12  2010/06/23 19:11:51  ssa1
#   Adding simple ribbon rendering and associated events for flexible fitting
#
#   Revision 1.11  2010/06/23 13:02:56  ssa1
#   Allowing users to reset a flexible fitting if need be.
#
#   Revision 1.10  2010/06/17 19:42:38  ssa1
#   Generic method for setting object specific coloring
#
#   Revision 1.9  2010/06/17 19:31:47  ssa1
#   Visually displaying flexible fitting clusters.
#
#   Revision 1.8  2010/05/27 05:19:31  ssa1
#   Moving all atoms when performing fitting instead of a single atom.
#
#   Revision 1.7  2010/05/26 21:53:21  ssa1
#   Adding in display styles for atom rendering.
#
#   Revision 1.6  2010/05/21 16:33:21  ssa1
#   Flexible fitting implemented in Gorgon
#
#   Revision 1.5  2010/05/21 16:11:45  ssa1
#   Flexible fitting implemented in Gorgon
#
#   Revision 1.4  2010/05/21 15:46:11  ssa1
#   Flexible fitting implemented in Gorgon
#
#   Revision 1.3  2010/05/21 15:45:16  ssa1
#   Flexible fitting implemented in Gorgon
#
#   Revision 1.2  2010/05/20 21:55:53  ssa1
#   Rigid body alignment based on largest flexible cluster
#
#   Revision 1.1  2010/05/20 19:15:15  ssa1
#   Flexible fitting interface.
#

import time
from math import sqrt, exp
from PyQt4 import QtCore, QtGui
from seq_model.Chain import Chain
from base_dock_widget import BaseDockWidget
from calpha_choose_chain_model import CAlphaChooseChainModel
from calpha_choose_chain_to_load_form import CAlphaChooseChainToLoadForm
from ui_dialog_calpha_flexible_fitting import Ui_DialogCAlphaFlexibleFitting
from libpyGORGON import Vector3DFloat, LinearSolver, SSECorrespondenceFinder, SSECorrespondenceFeature, FlexibleFittingEngine, MatrixFloat

class CAlphaFlexibleFittingForm(BaseDockWidget, Ui_DialogCAlphaFlexibleFitting):
        
    def __init__(self, main, viewer, parent=None):
        BaseDockWidget.__init__(self, 
                                main, 
                                "Fit to Volume", 
                                "Fit CAlpha Atoms to Volume", 
                                "perform_CAlphaFlexibleFitting", 
                                "actions-calpha-flexiblefitting", 
                                "actions-calpha", 
                                QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea | QtCore.Qt.BottomDockWidgetArea, 
                                QtCore.Qt.RightDockWidgetArea, 
                                parent)
        self.app = main

        self.calphaViewer   = viewer
        self.volumeViewer   = self.app.viewers["volume"]
        self.sseViewer      = self.app.viewers["sse"]
        self.skeletonViewer = self.app.viewers["skeleton"]

        self.linearSolver            = LinearSolver()
        self.sseCorrespondenceFinder = SSECorrespondenceFinder()
        self.flexibleFittingEngine   = FlexibleFittingEngine()

        self.sseViewer.renderer.setObjectSpecificColoring(True)

        self.alignToClusterStatus         = False
        self.previewCorrespondencesStatus = False

        self.createUI()

        # We need to load the correspondence
        self.correspondenceLoaded = False
        # pIndex, qIndex, isForward
        self.correspondence = {'pIndex':[], 'qIndex':[], 'isForward':[]}

        self.stage = "loadFiles"
        self.updateVisibility()

    def createUI(self):
        self.setupUi(self)    
       
        # Load Files
        self.connect(self.calphaViewer, QtCore.SIGNAL("modelLoaded()"), self.fileLoaded)
        self.connect(self.calphaViewer, QtCore.SIGNAL("modelUnloaded()"), self.fileLoaded)
        self.connect(self.sseViewer, QtCore.SIGNAL("modelLoaded()"), self.fileLoaded)
        self.connect(self.sseViewer, QtCore.SIGNAL("modelUnloaded()"), self.fileLoaded)
        self.connect(self.skeletonViewer, QtCore.SIGNAL("modelLoaded()"), self.fileLoaded)
        self.connect(self.skeletonViewer, QtCore.SIGNAL("modelUnloaded()"), self.fileLoaded)
        self.connect(self.pushButtonLoadFilesCAlphaAtoms, QtCore.SIGNAL("clicked (bool)"), self.loadCalphaAtoms)
        self.connect(self.pushButtonLoadFilesHelixAnnotations, QtCore.SIGNAL("clicked (bool)"), self.loadHelixAnnotations)
        self.connect(self.pushButtonLoadFilesSkeleton, QtCore.SIGNAL("clicked (bool)"), self.loadSkeleton)
        self.connect(self.pushButtonRigidDeformationPrincipalComponentTarget, QtCore.SIGNAL("clicked (bool)"), self.loadPrincipalComponentTargetChain)

        # For ground truth correspondence
        self.connect(self.pushButtonLoadCorrespondence, QtCore.SIGNAL("clicked (bool)"), self.loadCorrespondence)

        # Calculate Correspondence
        self.connect(self.pushButtonCorrespondenceOptionsCalculate, QtCore.SIGNAL("clicked (bool)"), self.calculateCorrespondences)

        # Alignment Selection
        self.connect(self.comboBoxAlignment, QtCore.SIGNAL("activated (int)"), self.updateCorrespondencesDisplay)

        # Align to Cluster
        self.connect(self.pushButtonCorrespondencesAlignToCluster, QtCore.SIGNAL("clicked (bool)"), self.alignToCluster)

        # Merge Clusters
        self.connect(self.pushButtonCorrespondencesMerge, QtCore.SIGNAL("clicked (bool)"), self.mergeClusters)

        # Manual Correspondence
        self.connect(self.pushButtonCorrespondencesManualCorrespondence, QtCore.SIGNAL("clicked (bool)"), self.manualCorrespondence)

        # Flip Correspondence
        self.connect(self.pushButtonCorrespondencesFlipCorrespondence, QtCore.SIGNAL("clicked (bool)"), self.flipCorrespondence)

        # Review Correspondence
        self.connect(self.pushButtonCorrespondencesPreview, QtCore.SIGNAL("clicked (bool)"), self.previewCorrespondences)
        self.connect(self.pushButtonCorrespondencesAccept, QtCore.SIGNAL("clicked (bool)"), self.acceptCorrespondences)
        self.connect(self.pushButtonCorrespondencesReject, QtCore.SIGNAL("clicked (bool)"), self.rejectCorrespondences)

        # Perform Deformation
        self.connect(self.pushButtonDeformationRigidDeform, QtCore.SIGNAL("clicked (bool)"), self.performRigidDeformation)
        self.connect(self.pushButtonDeformationFlexibleDeform, QtCore.SIGNAL("clicked (bool)"), self.performFlexibleDeformation)

        # Review Deformation
        self.connect(self.pushButtonDeformationAccept, QtCore.SIGNAL("clicked (bool)"), self.acceptDeformation)
        self.connect(self.pushButtonDeformationReject, QtCore.SIGNAL("clicked (bool)"), self.rejectDeformation)

    def updateVisibility(self):
        self.bringToFront()
        
        if self.stage == "loadFiles":
            self.pushButtonLoadFilesCAlphaAtoms.setEnabled(not self.calphaViewer.loaded)
            self.pushButtonLoadFilesHelixAnnotations.setEnabled(not self.sseViewer.loaded)
            self.pushButtonLoadFilesSkeleton.setEnabled(not self.skeletonViewer.loaded)
            # The correspondence does not come from any viewer.
            self.pushButtonLoadCorrespondence.setEnabled(not self.correspondenceLoaded)

            self.tabLoadFiles.setEnabled(True)

            if self.calphaViewer.loaded and self.sseViewer.loaded and self.skeletonViewer.loaded:
                self.tabCorrespondenceOptions.setEnabled(True)
            else:
                self.tabCorrespondenceOptions.setEnabled(False)

            self.tabCorrespondences.setEnabled(False)
            self.tabDeformation.setEnabled(False)
            self.tabDeformation.setEnabled(False)
            
            self.tabWidget.setCurrentWidget(self.tabLoadFiles)
        elif self.stage == "correspondenceOptions":
            self.tabLoadFiles.setEnabled(False)
            self.tabCorrespondenceOptions.setEnabled(True)
            self.tabCorrespondences.setEnabled(False)
            self.tabDeformationOptions.setEnabled(False)
            self.tabDeformation.setEnabled(False)
            
            self.tabWidget.setCurrentWidget(self.tabCorrespondenceOptions)
        elif self.stage == "correspondences":
            self.tabLoadFiles.setEnabled(False)
            self.tabCorrespondenceOptions.setEnabled(False)
            self.tabCorrespondences.setEnabled(True)
            self.tabDeformationOptions.setEnabled(False)
            self.tabDeformation.setEnabled(False)
            
            self.tabWidget.setCurrentWidget(self.tabCorrespondences)

            self.calphaViewer.setDisplayStyle(self.calphaViewer.DisplayStyleRibbon)
        elif self.stage == "deformationOptions":
            self.tabLoadFiles.setEnabled(False)
            self.tabCorrespondenceOptions.setEnabled(False)
            self.tabCorrespondences.setEnabled(False)
            self.tabDeformationOptions.setEnabled(True)
            self.tabDeformation.setEnabled(False)
            
            self.tabWidget.setCurrentWidget(self.tabDeformationOptions)
            
            self.calphaViewer.setDisplayStyle(self.calphaViewer.DisplayStyleRibbon)
            self.calphaViewer.setStrandVisibility(True)
            self.calphaViewer.setLoopVisibility(True)
        elif self.stage == "deformation":
            self.tabLoadFiles.setEnabled(False)
            self.tabCorrespondenceOptions.setEnabled(False)
            self.tabCorrespondences.setEnabled(False)
            self.tabDeformationOptions.setEnabled(False)
            self.tabDeformation.setEnabled(True)
            
            self.tabWidget.setCurrentWidget(self.tabDeformation)
            
            self.calphaViewer.setDisplayStyle(self.calphaViewer.DisplayStyleRibbon)
            self.calphaViewer.setStrandVisibility(True)
            self.calphaViewer.setLoopVisibility(True)

    def loadCorrespondence(self, temp):

        filename = unicode(QtGui.QFileDialog.getOpenFileName(self, self.tr("Open Data"), "", self.tr('Atom Positions (*.cor)')))
        extension = filename.split('.')[-1].lower()
        if extension == 'cor':
            for line in open(filename, 'U'):
                line = line.split(' ')
                if line[0] =='pIndex':
                    for i in line[1:]:
                        self.correspondence['pIndex'].append(int(i))
                elif line[0] == 'qIndex':
                    for i in line[1:]:
                        self.correspondence['qIndex'].append(int(i))
                elif line[0] == 'isForward':
                    for i in line[1:]:
                        self.correspondence['isForward'].append(bool(int(i)))
                else:
                    raise IOError, 'Invalid prefix in correspondence file...'
                    sys.exit(0)
        else:
            raise IOError, 'Wrong extention for correspondence file'
            sys.exit(0)

        self.correspondenceLoaded = True

        self.fileLoaded()

    def loadFixedCorrespondence(self, filename, dropHelicesNum):

        extension = filename.split('.')[-1].lower()
        if extension == 'cor':
            for line in open(filename, 'U'):
                line = line.split(' ')
                if line[0] =='pIndex':
                    for i in line[1:]:
                        self.correspondence['pIndex'].append(int(i))
                elif line[0] == 'qIndex':
                    for i in line[1:]:
                        self.correspondence['qIndex'].append(int(i))
                elif line[0] == 'isForward':
                    for i in line[1:]:
                        self.correspondence['isForward'].append(bool(int(i)))
                else:
                    raise IOError, 'Invalid prefix in correspondence file...'
                    sys.exit(0)
        else:
            raise IOError, 'Wrong extention for correspondence file'
            sys.exit(0)

        # Unifomrly drop the helices
        #print 'before:', self.correspondence['pIndex']
        n = len(self.correspondence['isForward'])
        k = dropHelicesNum
        # True means we choose to keep that match pair
        delta = n/k
        dropFlag = True
        mask = [True] * len(self.correspondence['isForward'])
        if delta < 2:
            dropFlag = False
            k = n-k

        remainOffset = n%(k+1)
        delta = n/(k+1)
        i = 0 if delta < 2 else delta

        while i < len(self.correspondence['isForward']) and k > 0:
            mask[i] = False
            i = i+delta+1 if remainOffset > 0 else i+delta
            remainOffset -= 1
            k -= 1

        if dropFlag:
            self.correspondence['isForward'] = [x for i, x in enumerate(self.correspondence['isForward']) if mask[i]]
            self.correspondence['pIndex'] = [x for i, x in enumerate(self.correspondence['pIndex']) if mask[i]]
            self.correspondence['qIndex'] = [x for i, x in enumerate(self.correspondence['qIndex']) if mask[i]]
        else:
            self.correspondence['isForward'] = [x for i, x in enumerate(self.correspondence['isForward']) if not mask[i]]
            self.correspondence['pIndex'] = [x for i, x in enumerate(self.correspondence['pIndex']) if not mask[i]]
            self.correspondence['qIndex'] = [x for i, x in enumerate(self.correspondence['qIndex']) if not mask[i]]

        #print 'mask:', mask
        #print 'after:', self.correspondence['pIndex']
        self.correspondenceLoaded = True
        self.fileLoaded()

    def loadCalphaAtoms(self, temp):
        self.app.actions.getAction("load_CAlpha").trigger()
        self.fileLoaded()
        
        self.chain = self.calphaViewer.loadedChains[0]

    def loadFixedCalphaAtoms(self, fileName, chainID):
        self.calphaViewer.loadFixedData(fileName, chainID)
        self.fileLoaded()

        self.chain = self.calphaViewer.loadedChains[0]

    def loadHelixAnnotations(self, temp):
        self.app.actions.getAction("load_SSE_Helix").trigger()
        self.fileLoaded()

    def loadFixedHelixAnnotations(self, fileName):
        #self.sseViewer.loadFixedHelixData(fileName)
        self.sseViewer.loadHelixDataFromFile(fileName)
        self.fileLoaded()

    def loadSkeleton(self, temp):
        self.app.actions.getAction("load_Skeleton").trigger()
        self.fileLoaded()

    def loadFixedSkeleton(self, fileName):
        self.skeletonViewer.loadDataFromFile(fileName)
        self.fileLoaded()

    def loadPrincipalComponentTargetChain(self, temp, filename=None, chainID=None):
        # Get Filename
        if not filename:
            filename = unicode(QtGui.QFileDialog.getOpenFileName(self, self.tr("Open Data"), "", self.tr('Atom Positions (*.pdb)\nFASTA (*.fas *.fa *.fasta)')))
        
        # Get Chain ID
        if not chainID:
            form = CAlphaChooseChainToLoadForm(filename)
            form.exec_()
            chainID = form.whichChainID

        # Load Chain
        if chainID != "ALL":
            self.principalComponentTargetChain = Chain.load(filename, None, chainID)
       
    def fileLoaded(self):
        if not self.calphaViewer.loaded or not self.sseViewer.loaded or not self.skeletonViewer.loaded or not self.correspondenceLoaded:
            self.stage = "loadFiles"
        else:
            self.stage = "correspondenceOptions"

        self.updateVisibility()

    def calculateCorrespondences(self, temp):
        # Create Correspondence Index Mapping
        self.correspondencePIndexToCalphaHelixIndex = {}
        
        # Load Calpha Helixes
        correspondenceIndex = 0
        for helixIndex, helix in self.chain.helices.items():
            # Get Atom Positions
            atomPositions = []
            for residueIndex in range(helix.startIndex, helix.stopIndex + 1):
                if residueIndex in self.chain.residueRange():
                    atomPositions.append(self.chain[residueIndex].getAtom('CA').getPosition())

            if atomPositions:
                # Calculate Best Fit Line
                startPosition = Vector3DFloat(0, 0, 0)
                endPosition   = Vector3DFloat(0, 0, 0)
                self.linearSolver.findBestFitLine(startPosition, endPosition, atomPositions)

                # Add to Correspondence Index Mapping 
                self.correspondencePIndexToCalphaHelixIndex[correspondenceIndex] = helixIndex
                correspondenceIndex += 1

                # Add to Engine
                self.flexibleFittingEngine.addCorrespondenceSourceFeature(startPosition, endPosition)

                # print helixIndex,' (start, end):',startPosition.x(),',',startPosition.y(),',',startPosition.z()
                # print '                ',endPosition.x(),',',endPosition.y(),',',endPosition.z()

        # Load Helix Annotations
        for helixIndex in range(self.sseViewer.renderer.getHelixCount()):
            # Get Volume Helix Line
            startPosition = self.sseViewer.renderer.getHelixCorner(helixIndex, 0)
            endPosition   = self.sseViewer.renderer.getHelixCorner(helixIndex, 1)

            # print helixIndex,' (start, end):',startPosition.x(),',',startPosition.y(),',',startPosition.z()
            # print '                ',endPosition.x(),',',endPosition.y(),',',endPosition.z()

            # Add to Engine
            self.flexibleFittingEngine.addCorrespondenceTargetFeature(startPosition, endPosition)

        # Initialize Correspondence Options
        self.flexibleFittingEngine.setCorrespondenceOptions(
            self.spinBoxMaximumAlignments.value(), 
            self.spinBoxHelixLengthDifference.value(), 
            self.spinBoxHelixCentroidDifference.value(), 
            self.spinBoxJointAngleDifference.value(), 
            self.spinBoxDihedralAngleDifference.value()
        )

        print "source helix", correspondenceIndex
        print "target helix", self.sseViewer.renderer.getHelixCount()

        # Set weighted svd weight range
        self.flexibleFittingEngine.setWSVDOptions(self.doubleSpinBoxBoudingBoxRatio.value())

        # Calculate Correspondences
        # Well, suppose we have done correspondences clear somewhere before
        # If loaded correspondence, then just use the loaded data : one alignment, one cluster and correspondence
        if self.correspondenceLoaded:
            for i in range(len(self.correspondence['pIndex'])):
                self.flexibleFittingEngine.addCorrespondencesByReadIn(self.correspondence['pIndex'][i], self.correspondence['qIndex'][i], self.correspondence['isForward'][i])
        else:
            self.flexibleFittingEngine.calculateCorrespondences()

        # print 'helix dir:', self.correspondence['isForward']

        # Update Correspondences Display
        self.updateCorrespondencesDisplay()

        # Update Visibility
        self.stage = "correspondences"
        self.updateVisibility()

    def updateCorrespondencesDisplay(self, temp = None):
        # Get Previous Alignment Index 
        previousAlignmentIndex = self.comboBoxAlignment.currentIndex() 
       
        # Calculate Renderer Calpha Helix Index Mapping
        rendererCalphaHelixIndexMapping = {}
        for rendererCalphaHelixIndex, calphaHelixIndex in enumerate(self.chain.helices.keys()):
            rendererCalphaHelixIndexMapping[calphaHelixIndex] = rendererCalphaHelixIndex

        # Set All Helixes to Unmatched Color
        unmatchedColor = QtGui.QColor()
        unmatchedColor.setRgbF(1.0, 1.0, 1.0)
        for calphaHelixIndex in self.chain.helices.keys():
            self.calphaViewer.renderer.setHelixColor(rendererCalphaHelixIndexMapping[calphaHelixIndex], unmatchedColor.redF(), unmatchedColor.greenF(), unmatchedColor.blueF())
        for volumeHelixIndex in range(self.sseViewer.renderer.getHelixCount()):
            self.sseViewer.renderer.setHelixColor(volumeHelixIndex, unmatchedColor.redF(), unmatchedColor.greenF(), unmatchedColor.blueF(), 1.0)
     
        # Clear Combo Box Alignment
        self.comboBoxAlignment.clear()
   
        # Clear Combo Box Cluster
        self.comboBoxCluster.clear()
        self.comboBoxFirstCluster.clear()
        self.comboBoxSecondCluster.clear()
   
        # Clear Table Widget Correspondences
        self.tableWidgetCorrespondences.clearContents()
        self.tableWidgetCorrespondences.setRowCount(0)
    
        # Set Header Resize Mode 
        self.tableWidgetCorrespondences.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
     
        # Get Alignment Count 
        alignmentCount = self.flexibleFittingEngine.getCorrespondenceAlignmentCount()
        if alignmentCount != 0:
            # Update Combo Box Alignment
            for alignmentIndex in range(alignmentCount):
                self.comboBoxAlignment.addItem(str(alignmentIndex))
            if previousAlignmentIndex != -1:
                self.comboBoxAlignment.setCurrentIndex(previousAlignmentIndex)
            
            # Get Selected Alignment Index
            selectedAlignmentIndex = int(self.comboBoxAlignment.currentText())

            # Update Table Widget Correspondences
            correspondenceIndex = 0
            clusterCount        = self.flexibleFittingEngine.getCorrespondenceClusterCount(selectedAlignmentIndex)
            for clusterIndex in range(clusterCount):
                # Update Combo Box Cluster
                self.comboBoxCluster.addItem(str(clusterIndex))
                self.comboBoxFirstCluster.addItem(str(clusterIndex))
                self.comboBoxSecondCluster.addItem(str(clusterIndex))

                # Cluster Color
                clusterColor = QtGui.QColor()
                clusterColor.setHsvF(float(clusterIndex) / float(clusterCount), 1.0, 1.0)
                clusterBrush = QtGui.QBrush(clusterColor)

                featureCount = self.flexibleFittingEngine.getCorrespondenceFeatureCount(selectedAlignmentIndex, clusterIndex)
                for featureIndex in range(featureCount):
                    # Get Correspondence
                    correspondence = self.flexibleFittingEngine.getCorrespondence(selectedAlignmentIndex, clusterIndex, featureIndex)
                   
                    # Get Helix Indices and Direction 
                    calphaHelixIndex = self.correspondencePIndexToCalphaHelixIndex[correspondence.getPIndex()]
                    volumeHelixIndex = correspondence.getQIndex()
                    isForward        = correspondence.isForward()

                    # Determine Direction
                    if isForward:
                        direction = "Forward"
                    else:
                        direction = "Backward"

                    # Update Helix Colors
                    tempColorScale = 3.0
                    self.calphaViewer.renderer.setHelixColor(rendererCalphaHelixIndexMapping[calphaHelixIndex], clusterColor.redF(), clusterColor.greenF(), clusterColor.blueF())
                    self.sseViewer.renderer.setHelixColor(volumeHelixIndex, tempColorScale*clusterColor.redF(), tempColorScale*clusterColor.greenF(), tempColorScale*clusterColor.blueF(), 1.0)

                    # print "helix color: ", calphaHelixIndex, volumeHelixIndex, ":", clusterColor.redF(), clusterColor.greenF(), clusterColor.blueF()

                    # Insert Table Row 
                    self.tableWidgetCorrespondences.insertRow(correspondenceIndex)
                   
                    # Create Table Items
                    clusterIndexItem     = QtGui.QTableWidgetItem(str(clusterIndex))
                    # Now the helix index display also starts from 0
                    calphaHelixIndexItem = QtGui.QTableWidgetItem(str(calphaHelixIndex-1))
                    volumeHelixIndexItem = QtGui.QTableWidgetItem(str(volumeHelixIndex))
                    directionItem        = QtGui.QTableWidgetItem(str(direction))

                    # Set Background Colors              
                    clusterIndexItem.setBackground(clusterBrush)
                    calphaHelixIndexItem.setBackground(clusterBrush)
                    volumeHelixIndexItem.setBackground(clusterBrush)
                    directionItem.setBackground(clusterBrush)
                    
                    # Set Table Items
                    self.tableWidgetCorrespondences.setItem(correspondenceIndex, 0, clusterIndexItem)
                    self.tableWidgetCorrespondences.setItem(correspondenceIndex, 1, calphaHelixIndexItem)
                    self.tableWidgetCorrespondences.setItem(correspondenceIndex, 2, volumeHelixIndexItem)
                    self.tableWidgetCorrespondences.setItem(correspondenceIndex, 3, directionItem)

                    correspondenceIndex += 1
    
        # Rerender
        self.calphaViewer.emitModelChanged()
        self.sseViewer.emitModelChanged()

    def manualCorrespondence(self, temp):
        # Get Manual Correspondence
        selectedAlignmentIndex     = int(self.comboBoxAlignment.currentText())
        selectedCalphaHelixIndices = self.calphaViewer.renderer.getSelectedHelixIndices()
        selectedVolumeHelixIndices = self.sseViewer.renderer.getSelectedHelixIndices()

        if len(selectedCalphaHelixIndices) == 1 and len(selectedVolumeHelixIndices) == 1:
            # Add Manual Correspondence
            self.flexibleFittingEngine.addCorrespondenceManual(selectedAlignmentIndex, selectedCalphaHelixIndices[0], selectedVolumeHelixIndices[0], False)

            # Update Correspondences Display
            self.updateCorrespondencesDisplay()

    def flipCorrespondence(self, temp):
        # Get Selected Alignment and Volume Helix Indices
        selectedAlignmentIndex     = int(self.comboBoxAlignment.currentText())
        selectedCalphaHelixIndices = self.calphaViewer.renderer.getSelectedHelixIndices()
        selectedVolumeHelixIndices = self.sseViewer.renderer.getSelectedHelixIndices()
        numSelectedHelixIndices    = len(selectedCalphaHelixIndices) + len(selectedVolumeHelixIndices)

        if len(selectedCalphaHelixIndices) == 1 and len(selectedVolumeHelixIndices) == 1:
            clusterCount = self.flexibleFittingEngine.getCorrespondenceClusterCount(selectedAlignmentIndex)
            for clusterIndex in range(clusterCount):
                featureCount = self.flexibleFittingEngine.getCorrespondenceFeatureCount(selectedAlignmentIndex, clusterIndex)
                for featureIndex in range(featureCount):
                    # Get Correspondence
                    correspondence = self.flexibleFittingEngine.getCorrespondence(selectedAlignmentIndex, clusterIndex, featureIndex)
                   
                    # Get Helix Indices and Direction
                    calphaHelixIndex = correspondence.getPIndex()
                    volumeHelixIndex = correspondence.getQIndex()
                    isForward        = correspondence.isForward()

                    # Flip Correspondence
                    if selectedCalphaHelixIndices[0] == calphaHelixIndex and selectedVolumeHelixIndices[0] == volumeHelixIndex:
                        self.flexibleFittingEngine.setCorrespondenceDirection(selectedAlignmentIndex, clusterIndex, featureIndex, not isForward)

            # Update Correspondences Display
            self.updateCorrespondencesDisplay()

    def alignToCluster(self, temp = None):
        # Disable Preview Correspondences
        if self.previewCorrespondencesStatus:
            self.previewCorrespondences()

        # Check Align To Cluster Status 
        if self.alignToClusterStatus:
            # Restore Chain
            self.restoreChain()
        else:
            # Backup Chain
            self.backupChain()

            # Get Selected Alignment Index
            selectedAlignmentIndex = int(self.comboBoxAlignment.currentText())
            
            # Get Selected Cluster Index
            selectedClusterIndex = int(self.comboBoxCluster.currentText())
            
            # Get Correspondence Transformation
            correspondenceTransformation = self.flexibleFittingEngine.getCorrespondenceClusterTransformation(selectedAlignmentIndex, selectedClusterIndex)

            # Show Align To Cluster 
            for residueIndex in self.chain.residueRange():
                for atomName in self.chain[residueIndex].getAtomNames():
                    self.chain[residueIndex].getAtom(atomName).transform(correspondenceTransformation)

        # Update Align To Cluster Status
        self.alignToClusterStatus = not self.alignToClusterStatus

        # Rerender
        self.calphaViewer.emitModelChanged()
        self.sseViewer.emitModelChanged()
   
    def mergeClusters(self, temp = None):
        # Get Selected Alignment Index
        selectedAlignmentIndex = int(self.comboBoxAlignment.currentText())
            
        # Get Selected Cluster Indices 
        selectedFirstClusterIndex  = int(self.comboBoxFirstCluster.currentText())
        selectedSecondClusterIndex = int(self.comboBoxSecondCluster.currentText())
 
        if selectedFirstClusterIndex != selectedSecondClusterIndex:
            # Merge Clusters 
            self.flexibleFittingEngine.mergeCorrespondenceClusters(selectedAlignmentIndex, selectedFirstClusterIndex, selectedSecondClusterIndex)
      
            # Update Correspondences Display
            self.updateCorrespondencesDisplay()
   
    def previewCorrespondences(self, temp = None):
        # Disable Align To Cluster
        if self.alignToClusterStatus:
            self.alignToCluster()

        # Check Preview Correspondences Status
        if self.previewCorrespondencesStatus:
            # Restore Chain
            self.restoreChain()
        else:
            # Backup Chain
            self.backupChain()

            # Calculate Missing Residue Indices 
            (missingResidueIndices, missingResidueIndexMapping) = self.calculateMissingResidueIndices()
        
            # Calculate Correspondence Transformations 
            (correspondences, correspondenceTransformations, correspondenceTargets) = self.calculateCorrespondenceTransformations()
            
            # Calculate Rigid Initialization Transformations
            rigidInitializationTransformations = self.calculateRigidInitializationTransformations(correspondences, correspondenceTransformations, correspondenceTargets, missingResidueIndexMapping)

            # Show Rigid Initialization
            for residueIndex in self.chain.residueRange():
                for atomName in self.chain[residueIndex].getAtomNames():
                    self.chain[residueIndex].getAtom(atomName).transform(rigidInitializationTransformations[residueIndex])

        # Update Preview Correspondences Status
        self.previewCorrespondencesStatus = not self.previewCorrespondencesStatus

        # Rerender
        self.calphaViewer.emitModelChanged()
        self.sseViewer.emitModelChanged()
        
    def acceptCorrespondences(self, temp):
        # Update Visibility
        self.stage = "deformationOptions"
        self.updateVisibility()

    def rejectCorrespondences(self, temp):
        # Disable Align To Cluster
        if self.alignToClusterStatus:
            self.alignToCluster()

        # Disable Preview Correspondeinces
        if self.previewCorrespondencesStatus:
            self.previewCorrespondences()
        
        # Reset Correspondences 
        self.flexibleFittingEngine.resetCorrespondence()

        # Update Correspondences Display
        self.updateCorrespondencesDisplay()
        
        # Update Visibility
        self.stage = "correspondenceOptions"
        self.updateVisibility()

    def performRigidDeformation(self, temp, selectedRigidFittingMethod=None):
        # Disable Align To Cluster
        if self.alignToClusterStatus:
            self.alignToCluster()

        # Disable Preview Correspondeinces
        if self.previewCorrespondencesStatus:
            self.previewCorrespondences()
       
        # Backup Chain
        self.backupChain()

        # Get Selected Alignment Index
        selectedAlignmentIndex = int(self.comboBoxAlignment.currentText())

        # Get Transformnation for Selected Rigid Fitting Method
        if not selectedRigidFittingMethod:
            selectedRigidFittingMethod = str(self.comboBoxDeformationRigidFittingMethod.currentText())
        if selectedRigidFittingMethod == "Alignment":
            # Calculate Alignment Transformation
            transformation = self.flexibleFittingEngine.getCorrespondenceAlignmentTransformation(selectedAlignmentIndex)
        elif selectedRigidFittingMethod == "Largest Cluster":
            # Get Largest Cluster Index
            largestClusterIndex = None
            for clusterIndex in range(self.flexibleFittingEngine.getCorrespondenceClusterCount(selectedAlignmentIndex)):
                clusterSize = self.flexibleFittingEngine.getCorrespondenceFeatureCount(selectedAlignmentIndex, clusterIndex)

                if largestClusterIndex == None or clusterSize > largestClusterSize:
                    largestClusterIndex = clusterIndex
                    largestClusterSize  = clusterSize

            # Calculate Largest Cluster Transformation
            transformation = self.flexibleFittingEngine.getCorrespondenceClusterTransformation(selectedAlignmentIndex, largestClusterIndex)
        elif selectedRigidFittingMethod == "Principal Components":
            # Calculate Principal Component Transformation
            transformation = self.calculatePrincipalComponentTransformation()

        # Transform Backbone and Side Chains
        for residueIndex in self.chain.residueRange():
            for atomName in self.chain[residueIndex].getAtomNames():
                self.chain[residueIndex].getAtom(atomName).transform(transformation)

        # Rerender
        self.calphaViewer.emitModelChanged()
        self.sseViewer.emitModelChanged()
        
        # Update Visibility
        self.stage = "deformation"
        self.updateVisibility()

    def performFlexibleDeformation(self, temp):
        # Disable Align To Cluster
        if self.alignToClusterStatus:
            self.alignToCluster()

        # Disable Preview Correspondeinces
        if self.previewCorrespondencesStatus:
            self.previewCorrespondences()
       
        # Backup Chain
        self.backupChain()
        
        # Calculate Missing Residue Indices
        (missingResidueIndices, missingResidueIndexMapping) = self.calculateMissingResidueIndices()

        # print "missingResidueIndices: ", missingResidueIndices

        # Calculate Missing Residue Bonds
        missingResidueBonds = []
        for currentResidueIndex, nextResidueIndex in self.ntuples(self.chain.residueRange(), 2):
            bond = (currentResidueIndex, nextResidueIndex)
            if bond not in self.chain.bonds and currentResidueIndex not in missingResidueIndices and nextResidueIndex not in missingResidueIndices:
                missingResidueBonds.append(bond)

        # Calculate Correspondence Transformations
        # self.calculateCorrespondenceTransformationsWithClusterWSVD()
        # (correspondences, correspondenceTransformations, correspondenceTargets) = self.calculateCorrespondenceTransformations()
        if self.checkBoxWSVD.isChecked():
            (correspondences, correspondenceTransformations, correspondenceTargets) = self.calculateCorrespondenceTransformationsWithClusterWSVD()
        else:
            (correspondences, correspondenceTransformations, correspondenceTargets) = self.calculateCorrespondenceTransformations()

        ### Global Fitting Options ###

        # Set Deformation Parameters
        self.flexibleFittingEngine.setDeformationNeighborhoodSizes(self.spinBoxLaplacianNeighborhoodSize.value(), self.spinBoxTransformationNeighborhoodSize.value(), 3)
        
        # Add Deformation Original Vertices
        for residueIndex in self.chain.residueRange():
            if residueIndex not in missingResidueIndices: 
                self.flexibleFittingEngine.addDeformationOriginalVertex(self.chain[residueIndex].getAtom('CA').getPosition())
                # print self.chain[residueIndex].getAtom('CA').getPosition().x()

        # Add Deformation Edges and Calculate Neighbor Vertex Indices Cache
        numResidues      = len(self.chain) - len(missingResidueIndices)
        deformationEdges = set()
        for bond in self.chain.bonds:
            deformationEdges.add((missingResidueIndexMapping[bond[0]], missingResidueIndexMapping[bond[1]]))
            #if bond[0] == 130 or bond[1] == 130 or bond[0] == 129 or bond[1] == 129 or bond[0] == 131 or bond[1] == 131:
                #tbs = self.chain[bond[0]].getAtom('CA').getPosition()
                #tbe = self.chain[bond[1]].getAtom('CA').getPosition()
                #print "helix", tbs.x(),",",tbs.y(),",",tbs.z(),",",tbe.x(),",",tbe.y(),",",tbe.z()
            #if bond[0] == 252 or bond[1] == 252 or bond[0] == 251 or bond[1] == 251 or bond[0] == 253 or bond[1] == 253:
                #tbs = self.chain[bond[0]].getAtom('CA').getPosition()
                #tbe = self.chain[bond[1]].getAtom('CA').getPosition()
                #print "loop: ", tbs.x(),",",tbs.y(),",",tbs.z(),",",tbe.x(),",",tbe.y(),",",tbe.z()
            #if bond[0] == 172 or bond[1] == 172 or bond[0] == 171 or bond[1] == 171 or bond[0] == 173 or bond[1] == 173 or bond[0] == 142 or bond[1] == 142 or bond[0] == 140 or bond[1] == 140:
                #tbs = self.chain[bond[0]].getAtom('CA').getPosition()
                #tbe = self.chain[bond[1]].getAtom('CA').getPosition()
                #print "sheet: ", tbs.x(),",",tbs.y(),",",tbs.z(),",",tbe.x(),",",tbe.y(),",",tbe.z()

        for bond in missingResidueBonds:
            deformationEdges.add((missingResidueIndexMapping[bond[0]], missingResidueIndexMapping[bond[1]]))

        for sheet in self.chain.sheets.values():
            for bond in sheet.bonds:
                deformationEdges.add((missingResidueIndexMapping[bond[0]], missingResidueIndexMapping[bond[1]]))
                #tbs = self.chain[bond[0]].getAtom('CA').getPosition()
                #tbe = self.chain[bond[1]].getAtom('CA').getPosition()
                #print "sheet: ", tbs.x(),",",tbs.y(),",",tbs.z(),",",tbe.x(),",",tbe.y(),",",tbe.z()

        for deformationEdge in deformationEdges:
            self.flexibleFittingEngine.addDeformationEdge(deformationEdge[0], deformationEdge[1])

        self.flexibleFittingEngine.calculateNeighborVertexIndicesCaches(numResidues)

        ### Start Helix Fitting ###

        # Set Deformation Vertices
        self.flexibleFittingEngine.setDeformationVerticesToOriginalVertices()

        # Only to test fitting weight and distortion weigh effect
        #self.flexibleFittingEngine.setDeformationWeights(5.0, 1.0)

        start = time.clock()

        # Calculate and Add Rigid Initialization Transformations
        self.flexibleFittingEngine.setDeformationRigidInitialization(self.checkBoxRigidInitialization.isChecked())
        if self.checkBoxRigidInitialization.isChecked():
            rigidInitializationTransformations = self.calculateRigidInitializationTransformations(correspondences, correspondenceTransformations, correspondenceTargets, missingResidueIndexMapping)
            
            for residueIndex in self.chain.residueRange():
                if residueIndex not in missingResidueIndices:
                    self.flexibleFittingEngine.addDeformationRigidInitializationTransformation(rigidInitializationTransformations[residueIndex])

        # Add Deformation Handles
        for calphaHelixIndex in sorted(correspondenceTargets.keys()):
            #print calphaHelixIndex
            for residueIndex in sorted(correspondenceTargets[calphaHelixIndex].keys()):
                self.flexibleFittingEngine.addDeformationHandle(missingResidueIndexMapping[residueIndex], correspondenceTargets[calphaHelixIndex][residueIndex])
                #tv = self.chain[residueIndex].getAtom('CA').getPosition()
                #print tv.x(),",",tv.y(),",",tv.z()

        # Perform Deformation
        self.flexibleFittingEngine.deformLaplacian()
        #self.flexibleFittingEngine.deformWSVD()

        elapsed = (time.clock() - start)
        print "helix-guided time: ", elapsed

        ### End Helix Fitting ###

        ### Start Skeleton Fitting ###
        
        # Set Deformation Parameters
        self.flexibleFittingEngine.setDeformationRigidInitialization(False)
        self.flexibleFittingEngine.setDeformationWeights(self.doubleSpinBoxSkeletonFittingWeight.value(), self.doubleSpinBoxSkeletonDistortionWeight.value())
  
        # Iterative Skeleton Fitting
        skeletonMesh                       = self.skeletonViewer.renderer.getMesh()
        originalSkeletonSurfaceVertices    = skeletonMesh.getSurfaceVertices(4)
        originalSkeletonNonSurfaceVertices = skeletonMesh.getNonSurfaceVertices()
        skeletonFittingIterations          = self.spinBoxSkeletonFittingIterations.value()
        skeletonFittingDistanceThreshold   = self.spinBoxSkeletonFittingDistanceThreshold.value()

        realIterationRound = 0;
        start = time.clock()

        # skeletonFittingIterations serves as the most iteration round if "auto stop" box is checked
        for skeletonFittingIteration in range(skeletonFittingIterations):

            realIterationRound +=1

            # Get Deformed Vertices
            deformedCalphaPositions = self.flexibleFittingEngine.getDeformedVertices()

            # Set Deformation Vertices
            self.flexibleFittingEngine.setDeformationVerticesToDeformedVertices()
        
            # Reset Deformation Handles
            self.flexibleFittingEngine.resetDeformationHandles()
 
            # Copy Original Skeleton Vertices
            skeletonSurfaceVertices    = list(originalSkeletonSurfaceVertices)
            skeletonNonSurfaceVertices = list(originalSkeletonNonSurfaceVertices)
  
            # Add Skeleton Deformation Handles
            for residueIndex in self.chain.residueRange():
                if residueIndex not in missingResidueIndices:
                    isHelix = False
                    isSheet = False
                    for helixIndex, helix in self.chain.helices.items():
                        if helixIndex in correspondenceTargets and helix.startIndex <= residueIndex and residueIndex <= helix.stopIndex:
                            isHelix                  = True
                            correspondenceHelixIndex = helixIndex
                    for sheet in self.chain.sheets.values():
                        for strand in sheet.strandList.values():
                            if strand.startIndex <= residueIndex and residueIndex <= strand.stopIndex:
                                isSheet = True

                    calphaAtomPosition = deformedCalphaPositions[missingResidueIndexMapping[residueIndex]]
                   
                    if isHelix:
                        closestSkeletonVertex = correspondenceTargets[correspondenceHelixIndex][residueIndex]
                    elif isSheet:
                        closestSkeletonVertex = self.getClosestVertex(calphaAtomPosition, skeletonSurfaceVertices, skeletonFittingDistanceThreshold, False)
                        if closestSkeletonVertex == None:
                            closestSkeletonVertex = self.getClosestVertex(calphaAtomPosition, skeletonNonSurfaceVertices, skeletonFittingDistanceThreshold, False)
                    else:
                        closestSkeletonVertex = self.getClosestVertex(calphaAtomPosition, skeletonNonSurfaceVertices, skeletonFittingDistanceThreshold, False)

                    if closestSkeletonVertex != None:
                        self.flexibleFittingEngine.addDeformationHandle(missingResidueIndexMapping[residueIndex], closestSkeletonVertex)
                        
                        #if skeletonFittingIteration == (skeletonFittingIterations - 1):
                        #    self.calphaViewer.renderer.addSkeletonFittingVector(self.chain[residueIndex].getAtom('CA').getPosition(), closestSkeletonVertex)
            
            # Perform Deformation
            self.flexibleFittingEngine.deformLaplacian()

            # If RMSD difference between deformed and pre-deformed model is below the threshold, stop
            if self.checkBoxAutoStop.isChecked() and self.flexibleFittingEngine.getShapeDistance() < self.doubleSpinBoxStopThreshold.value():
                break

        elapsed = (time.clock() - start)
        print "skeleton-guided time: ", elapsed
        print "skeleton iteration round: ", realIterationRound

        ### End Skeleton Fitting ###

        '''
        ### Start Final Helix Fitting --- Enforce helix overlap --- deprecated ###
        ### Right way to do: maybe we need adaptive distortion term for each atom ###
        ### High fitting for helix atoms and high distortion for others ###
        ### Maybe we only deform helix atoms ###

        self.flexibleFittingEngine.setDeformationRigidInitialization(False)
        currentFittingWeight = 0.9
        tempIteration = 2
        weightDelta = (1.0 - currentFittingWeight) / float(tempIteration)
        weightEpsilon = 0.01
        for i in range(tempIteration):
            # Set Fitting and Distortion weight
            self.flexibleFittingEngine.setDeformationWeights(currentFittingWeight, 1.0 - currentFittingWeight + weightEpsilon)

            # Set Deformation Vertices
            self.flexibleFittingEngine.setDeformationVerticesToDeformedVertices()

            # Clear Deformation Handles
            self.flexibleFittingEngine.resetDeformationHandles()
            # Set Deformation Handles
            for calphaHelixIndex in sorted(correspondenceTargets.keys()):
                for residueIndex in sorted(correspondenceTargets[calphaHelixIndex].keys()):
                    self.flexibleFittingEngine.addDeformationHandle(missingResidueIndexMapping[residueIndex], correspondenceTargets[calphaHelixIndex][residueIndex])

            # Perform Deformation
            self.flexibleFittingEngine.deformLaplacian()
            # Update Deformation weight
            currentFittingWeight += weightDelta

        ### End Final Helix Fitting ###
        '''

        start = time.clock()

        # Calculate Deformation Transformations
        self.flexibleFittingEngine.calculateDeformationTransformations(numResidues)

        # Transform Backbone and Side Chains
        for residueIndex in self.chain.residueRange():
            deformationTransformation = self.flexibleFittingEngine.getDeformationTransformation(missingResidueIndexMapping[residueIndex])
            for atomName in self.chain[residueIndex].getAtomNames():
                self.chain[residueIndex].getAtom(atomName).transform(deformationTransformation)

        elapsed = (time.clock() - start)
        print "side-chain time: ", elapsed

        # Rerender
        self.calphaViewer.emitModelChanged()
        self.sseViewer.emitModelChanged()
        
        # Update Visibility
        self.stage = "deformation"
        self.updateVisibility()

    def acceptDeformation(self, temp):
        # Clear Calpha Helix Colors
        self.calphaViewer.renderer.clearHelixColors()

        unmatchedColor = QtGui.QColor()
        unmatchedColor.setRgbF(0.0, 1.0, 0.0)
        for volumeHelixIndex in range(self.sseViewer.renderer.getHelixCount()):
            self.sseViewer.renderer.setHelixColor(volumeHelixIndex, unmatchedColor.redF(), unmatchedColor.greenF(), unmatchedColor.blueF(), 1.0)

        # Rerender
        self.calphaViewer.emitModelChanged()
        self.sseViewer.emitModelChanged()

        self.close()

    def rejectDeformation(self, temp):
        # Restore Chain
        self.restoreChain()

        # Reset Flexible Fitting Engine
        self.flexibleFittingEngine.resetDeformation()

        # Rerender
        self.calphaViewer.emitModelChanged()
        self.sseViewer.emitModelChanged()
        
        # Update Visibility
        self.stage = "deformationOptions"
        self.updateVisibility()

    def calculateMissingResidueIndices(self):
        missingResidueIndices      = {}
        missingResidueIndexMapping = {}

        mappedResidueIndex = 0
        for originalResidueIndex in self.chain.residueRange():
            calphaAtom = self.chain[originalResidueIndex].getAtom('CA')
            if calphaAtom:
                missingResidueIndexMapping[originalResidueIndex] = mappedResidueIndex
                mappedResidueIndex += 1
            else:
                missingResidueIndices[originalResidueIndex] = True

        for missingResidueIndex in missingResidueIndices.keys():
            for direction in (-1, 1):
                residueIndexIterator = missingResidueIndex
                while (residueIndexIterator >= self.chain.getFirstResidueIndex()):
                    residueIndexIterator += direction
                    if residueIndexIterator in missingResidueIndexMapping:
                        missingResidueIndexMapping[missingResidueIndex] = missingResidueIndexMapping[residueIndexIterator]
                        break

                if missingResidueIndex in missingResidueIndexMapping:
                    break

        return (missingResidueIndices, missingResidueIndexMapping)

    def calculatePrincipalComponentTransformation(self):
        # Options
        numAxes = 3

        # Get Source and Target Calpha Atoms
        sourceChain = self.chain
        targetChain = self.principalComponentTargetChain 

        # Check Target Chain
        if targetChain == None:
            return

        # Calculate Source Centroid and Number of Residues
        sourceCentroid    = Vector3DFloat(0.0, 0.0, 0.0)
        numSourceResidues = 0
        for residueIndex in sourceChain.residueRange():
            sourceCalphaAtom = sourceChain[residueIndex].getAtom('CA')
            if sourceCalphaAtom:
                sourceCentroid += sourceCalphaAtom.getPosition()
                numSourceResidues += 1
        sourceCentroid = sourceCentroid * (1.0 / numSourceResidues)

        # Construct Source Position Variances
        sourcePositionVariances = MatrixFloat(3, numSourceResidues)
        for residueNum, residueIndex in enumerate(sourceChain.residueRange()):
            sourceCalphaAtom = sourceChain[residueIndex].getAtom('CA')
            if sourceCalphaAtom:
                positionVariance = sourceCalphaAtom.getPosition() - sourceCentroid
                sourcePositionVariances.setValue(positionVariance.x(), 0, residueNum) 
                sourcePositionVariances.setValue(positionVariance.y(), 1, residueNum) 
                sourcePositionVariances.setValue(positionVariance.z(), 2, residueNum)

        # Perform Source Principal Component Analysis
        sourceAxes             = MatrixFloat(3, 3)
        sourceAxesRankings     = MatrixFloat(3, numSourceResidues)
        sourceAxesCoefficients = MatrixFloat(numSourceResidues, numSourceResidues)
        sourcePositionVariances.singularValueDecomposition(sourceAxes, sourceAxesRankings, sourceAxesCoefficients)

        # Compute Third Souce Axis
        sourceFirstAxis  = Vector3DFloat(sourceAxes.getValue(0, 0), sourceAxes.getValue(1, 0), sourceAxes.getValue(2, 0))
        sourceSecondAxis = Vector3DFloat(sourceAxes.getValue(0, 1), sourceAxes.getValue(1, 1), sourceAxes.getValue(2, 1))
        sourceThirdAxis  = sourceFirstAxis ^ sourceSecondAxis
        sourceThirdAxis.normalize()
        sourceAxes.setValue(sourceThirdAxis.x(), 0, 2)
        sourceAxes.setValue(sourceThirdAxis.y(), 1, 2)
        sourceAxes.setValue(sourceThirdAxis.z(), 2, 2)

        # Calculate Target Centroid and Number of Residues
        targetCentroid    = Vector3DFloat(0.0, 0.0, 0.0)
        numTargetResidues = 0
        for residueIndex in targetChain.residueRange():
            targetCalphaAtom = targetChain[residueIndex].getAtom('CA')
            if targetCalphaAtom:
                targetCentroid += targetCalphaAtom.getPosition()
                numTargetResidues += 1
        targetCentroid = targetCentroid * (1.0 / numTargetResidues)
        
        # Construct Target Position Variances
        targetPositionVariances = MatrixFloat(3, numTargetResidues)
        for residueNum, residueIndex in enumerate(targetChain.residueRange()):
            targetCalphaAtom = targetChain[residueIndex].getAtom('CA')
            if targetCalphaAtom:
                positionVariance = targetCalphaAtom.getPosition() - targetCentroid
                targetPositionVariances.setValue(positionVariance.x(), 0, residueNum) 
                targetPositionVariances.setValue(positionVariance.y(), 1, residueNum) 
                targetPositionVariances.setValue(positionVariance.z(), 2, residueNum)

        # Perform Target Principal Component Analysis
        targetAxes             = MatrixFloat(3, 3)
        targetAxesRankings     = MatrixFloat(3, numTargetResidues)
        targetAxesCoefficients = MatrixFloat(numTargetResidues, numTargetResidues)
        targetPositionVariances.singularValueDecomposition(targetAxes, targetAxesRankings, targetAxesCoefficients)

        # Compute Third Target Axis
        targetFirstAxis  = Vector3DFloat(targetAxes.getValue(0, 0), targetAxes.getValue(1, 0), targetAxes.getValue(2, 0))
        targetSecondAxis = Vector3DFloat(targetAxes.getValue(0, 1), targetAxes.getValue(1, 1), targetAxes.getValue(2, 1))
        targetThirdAxis  = targetFirstAxis ^ targetSecondAxis
        targetThirdAxis.normalize()
        targetAxes.setValue(targetThirdAxis.x(), 0, 2)
        targetAxes.setValue(targetThirdAxis.y(), 1, 2)
        targetAxes.setValue(targetThirdAxis.z(), 2, 2)

        # Translation To Origin
        translationToOrigin = MatrixFloat(4, 4)
        translationToOrigin.setValue(1.0                , 0, 0)
        translationToOrigin.setValue(1.0                , 1, 1)
        translationToOrigin.setValue(1.0                , 2, 2)
        translationToOrigin.setValue(1.0                , 3, 3)
        translationToOrigin.setValue(-sourceCentroid.x(), 0, 3)
        translationToOrigin.setValue(-sourceCentroid.y(), 1, 3)
        translationToOrigin.setValue(-sourceCentroid.z(), 2, 3)

        # Translation From Origin
        translationFromOrigin = MatrixFloat(4, 4)
        translationFromOrigin.setValue(1.0               , 0, 0)
        translationFromOrigin.setValue(1.0               , 1, 1)
        translationFromOrigin.setValue(1.0               , 2, 2)
        translationFromOrigin.setValue(1.0               , 3, 3)
        translationFromOrigin.setValue(targetCentroid.x(), 0, 3)
        translationFromOrigin.setValue(targetCentroid.y(), 1, 3)
        translationFromOrigin.setValue(targetCentroid.z(), 2, 3)

        # Get Selected Orientation
        selectedOrientation = self.comboBoxRigidDeformationPrincipalComponentOrientation.currentText()
        orientations        = [[1, 1, 1], [1, -1, -1], [-1, -1, 1], [-1, 1, -1]]
        if selectedOrientation == "1":
            orientations = [orientations[0]]
        elif selectedOrientation == "2":
            orientations = [orientations[1]]
        elif selectedOrientation == "3":
            orientations = [orientations[2]]
        elif selectedOrientation == "4":
            orientations = [orientations[3]]

        # Rotation
        bestRotation = None
        for orientation in orientations:
            # Construct Transposed Oriented Source Axes
            orientedSourceAxes = MatrixFloat(3, 3)
            for axisIndex in range(numAxes):
                orientedSourceAxes.setValue(orientation[axisIndex] * sourceAxes.getValue(0, axisIndex), 0, axisIndex)
                orientedSourceAxes.setValue(orientation[axisIndex] * sourceAxes.getValue(1, axisIndex), 1, axisIndex)
                orientedSourceAxes.setValue(orientation[axisIndex] * sourceAxes.getValue(2, axisIndex), 2, axisIndex)
            
            # Calculate Rotation and Angle
            rotation = targetAxes * orientedSourceAxes.transpose()

            # Convert Rotation to Homogeneous Coordinates
            homogeneousRotation = MatrixFloat(4 ,4)
            for axisIndex in range(numAxes):
                homogeneousRotation.setValue(rotation.getValue(0, axisIndex), 0, axisIndex)
                homogeneousRotation.setValue(rotation.getValue(1, axisIndex), 1, axisIndex)
                homogeneousRotation.setValue(rotation.getValue(2, axisIndex), 2, axisIndex)
            homogeneousRotation.setValue(1.0, 3, 3)

            # Get Transformed Source Positions
            transformedSourcePositions = []
            for residueIndex in sourceChain.residueRange():
                sourceCalphaAtom = sourceChain[residueIndex].getAtom('CA')
                if sourceCalphaAtom:
                    transformedSourcePositions.append(sourceCalphaAtom.getPosition().Transform(translationFromOrigin * homogeneousRotation * translationToOrigin))

            # Get Target Positions
            targetPositions = []
            for residueIndex in targetChain.residueRange():
                targetCalphaAtom = targetChain[residueIndex].getAtom('CA')
                if targetCalphaAtom:
                    targetPositions.append(targetCalphaAtom.getPosition())

            # Calculate RMSD
            rmsd = self.calculateRMSD(transformedSourcePositions, targetPositions, None, False)
            
            # Check Against Best RMSD
            if bestRotation == None or rmsd < bestRMSD:
                bestRotation = homogeneousRotation
                bestRMSD     = rmsd 
            
        # Calcuate Transformation
        principalComponentTransformation = translationFromOrigin * bestRotation * translationToOrigin

        return principalComponentTransformation

    def calculateCorrespondenceTransformations(self):
        correspondences               = {}
        correspondenceTransformations = {}
        correspondenceTargets         = {}
        selectedAlignmentIndex = int(self.comboBoxAlignment.currentText())
        clusterCount = self.flexibleFittingEngine.getCorrespondenceClusterCount(selectedAlignmentIndex)
        for clusterIndex in range(clusterCount):
            featureCount = self.flexibleFittingEngine.getCorrespondenceFeatureCount(selectedAlignmentIndex, clusterIndex)
            for featureIndex in range(featureCount):
                # Get Correspondence
                correspondence = self.flexibleFittingEngine.getCorrespondence(selectedAlignmentIndex, clusterIndex, featureIndex)
                
                # Get Indices
                calphaHelixIndex = self.correspondencePIndexToCalphaHelixIndex[correspondence.getPIndex()]
                volumeHelixIndex = correspondence.getQIndex()

                # Get Transformation
                correspondenceTransformation = self.flexibleFittingEngine.getCorrespondenceFeatureTransformation(selectedAlignmentIndex, clusterIndex, featureIndex)

                # Cache Correpondence
                correspondences[calphaHelixIndex] = volumeHelixIndex

                # Cache Correspondence Transformation
                correspondenceTransformations[calphaHelixIndex] = correspondenceTransformation
               
                # Cache Correspondence Targets
                helix = self.chain.helices[calphaHelixIndex]
                for residueIndex in range(helix.startIndex, helix.stopIndex + 1):
                    if residueIndex in self.chain.residueRange():
                        if calphaHelixIndex not in correspondenceTargets:
                            correspondenceTargets[calphaHelixIndex] = {}

                        correspondenceTargets[calphaHelixIndex][residueIndex] = self.chain[residueIndex].getAtom('CA').getPosition().Transform(correspondenceTransformation)

        return (correspondences, correspondenceTransformations, correspondenceTargets)

    def calculateCorrespondenceTransformationsWithClusterWSVD(self):
        correspondences               = {}
        correspondenceTransformations = {}
        correspondenceTargets         = {}
        selectedAlignmentIndex = int(self.comboBoxAlignment.currentText())
        clusterCount = self.flexibleFittingEngine.getCorrespondenceClusterCount(selectedAlignmentIndex)

        # Precompute the weight list for each correspondence
        self.flexibleFittingEngine.computeCorrespondencePairWiseWeight(selectedAlignmentIndex)

        for clusterIndex in range(clusterCount):
            featureCount = self.flexibleFittingEngine.getCorrespondenceFeatureCount(selectedAlignmentIndex, clusterIndex)
            for featureIndex in range(featureCount):
                # Get Correspondence
                correspondence = self.flexibleFittingEngine.getCorrespondence(selectedAlignmentIndex, clusterIndex, featureIndex)

                # Get Indices
                calphaHelixIndex = self.correspondencePIndexToCalphaHelixIndex[correspondence.getPIndex()]
                volumeHelixIndex = correspondence.getQIndex()

                # Get Transformation
                correspondenceTransformation = self.flexibleFittingEngine.getCorrespondenceFeatureTransformationWSVD(selectedAlignmentIndex, clusterIndex, featureIndex)

                # Cache Correpondence
                correspondences[calphaHelixIndex] = volumeHelixIndex

                # Cache Correspondence Transformation
                correspondenceTransformations[calphaHelixIndex] = correspondenceTransformation

                # Cache Correspondence Targets
                helix = self.chain.helices[calphaHelixIndex]
                for residueIndex in range(helix.startIndex, helix.stopIndex + 1):
                    if residueIndex in self.chain.residueRange():
                        if calphaHelixIndex not in correspondenceTargets:
                            correspondenceTargets[calphaHelixIndex] = {}

                        correspondenceTargets[calphaHelixIndex][residueIndex] = self.chain[residueIndex].getAtom('CA').getPosition().Transform(correspondenceTransformation)

        return (correspondences, correspondenceTransformations, correspondenceTargets)

    def calculateRigidInitializationTransformations(self, correspondences, correspondenceTransformations, correspondenceTargets, missingResidueIndexMapping):
        rigidInitializationTransformations = {}
        
        correspondedCalphaHelixIndices = sorted(correspondences.keys())
        
        # Loops
        for currentCalphaHelixIndex, nextCalphaHelixIndex in self.ntuples(correspondedCalphaHelixIndices, 2):
            # Get Helices
            currentCalphaHelix = self.chain.helices[currentCalphaHelixIndex]
            nextCalphaHelix    = self.chain.helices[nextCalphaHelixIndex]

            # Get Source/Target Positions
            sourceCalphaHelixPositions = []
            targetCalphaHelixPositions = []
            for currentResidueIndex in range(currentCalphaHelix.startIndex, currentCalphaHelix.stopIndex + 1):
                if currentResidueIndex in self.chain.residueRange():
                    sourceCalphaHelixPositions.append(self.chain[currentResidueIndex].getAtom('CA').getPosition())
                    targetCalphaHelixPositions.append(correspondenceTargets[currentCalphaHelixIndex][currentResidueIndex])
            for nextResidueIndex in range(nextCalphaHelix.startIndex, nextCalphaHelix.stopIndex + 1):
                if nextResidueIndex in self.chain.residueRange():
                    sourceCalphaHelixPositions.append(self.chain[nextResidueIndex].getAtom('CA').getPosition())
                    targetCalphaHelixPositions.append(correspondenceTargets[nextCalphaHelixIndex][nextResidueIndex])

            # Calculate Rigid Initialization Transformation
            rigidInitializationTransformation = self.linearSolver.findRotationTranslation(sourceCalphaHelixPositions, targetCalphaHelixPositions)

            # Store Rigid Initialization Transformation
            if currentCalphaHelixIndex == correspondedCalphaHelixIndices[0]:
                startResidueIndex = self.chain.getFirstResidueIndex()
            else:
                startResidueIndex = currentCalphaHelix.startIndex

            if nextCalphaHelixIndex == correspondedCalphaHelixIndices[-1]:
                stopResidueIndex = self.chain.getLastResidueIndex() + 2
            else:
                stopResidueIndex = nextCalphaHelix.startIndex

            for residueIndex in range(startResidueIndex, stopResidueIndex + 1):
                if residueIndex in self.chain.residueRange():
                    rigidInitializationTransformations[residueIndex] = rigidInitializationTransformation

        # Helices
        for calphaHelixIndex in correspondedCalphaHelixIndices:
            # Get Helix
            calphaHelix = self.chain.helices[calphaHelixIndex]

            # Get Helix Transformation
            helixTransformation = correspondenceTransformations[calphaHelixIndex]

            # Store Rigid Initialization Transformation
            for residueIndex in range(calphaHelix.startIndex, calphaHelix.stopIndex + 1):
                if residueIndex in self.chain.residueRange():
                    rigidInitializationTransformations[residueIndex] = helixTransformation

        # Sheets
        for sheetNo, sheet in sorted(self.chain.sheets.items()):
            # Get Connected Helices
            connectedCalphaHelixIndices = set()
            for strandNo, strand in sheet.strandList.items():
                for currentCalphaHelixIndex, nextCalphaHelixIndex in self.ntuples(correspondedCalphaHelixIndices, 2):
                    currentCalphaHelix = self.chain.helices[currentCalphaHelixIndex]
                    nextCalphaHelix    = self.chain.helices[nextCalphaHelixIndex]

                    if currentCalphaHelixIndex == correspondedCalphaHelixIndices[0] and strand.startIndex >= self.chain.getFirstResidueIndex() and strand.stopIndex < nextCalphaHelix.stopIndex:
                        connectedCalphaHelixIndices.add(currentCalphaHelixIndex)
                        connectedCalphaHelixIndices.add(nextCalphaHelixIndex)
                    
                    if nextCalphaHelixIndex == correspondedCalphaHelixIndices[-1] and strand.startIndex > currentCalphaHelix.startIndex and strand.stopIndex <= self.chain.getLastResidueIndex():
                        connectedCalphaHelixIndices.add(currentCalphaHelixIndex)
                        connectedCalphaHelixIndices.add(nextCalphaHelixIndex)

                    if strand.startIndex > currentCalphaHelix.stopIndex and strand.stopIndex < nextCalphaHelix.startIndex:
                        connectedCalphaHelixIndices.add(currentCalphaHelixIndex)
                        connectedCalphaHelixIndices.add(nextCalphaHelixIndex)
            if not connectedCalphaHelixIndices:
                continue

            # Get Source/Target Positions
            sourceCalphaHelixPositions = []
            targetCalphaHelixPositions = []
            for calphaHelixIndex in connectedCalphaHelixIndices:
                calphaHelix = self.chain.helices[calphaHelixIndex]

                for residueIndex in range(calphaHelix.startIndex, calphaHelix.stopIndex + 1):
                    if residueIndex in self.chain.residueRange():
                        sourceCalphaHelixPositions.append(self.chain[residueIndex].getAtom('CA').getPosition())
                        targetCalphaHelixPositions.append(correspondenceTargets[calphaHelixIndex][residueIndex])

            # Calculate Rigid Initialization Transformation
            rigidInitializationTransformation = self.linearSolver.findRotationTranslation(sourceCalphaHelixPositions, targetCalphaHelixPositions)

            # Store Rigid Initialization Transformation
            for strand in sheet.strandList.values():
                for residueIndex in range(strand.startIndex, strand.stopIndex + 1):
                    if residueIndex in self.chain.residueRange():
                        rigidInitializationTransformations[residueIndex] = rigidInitializationTransformation

        return rigidInitializationTransformations

    def backupChain(self):
        self.backupChainAtomPositions = {}
        for residueIndex in self.chain.residueRange():
            if residueIndex not in self.backupChainAtomPositions:
                self.backupChainAtomPositions[residueIndex] = {}

            for atomName in self.chain[residueIndex].getAtomNames():
                self.backupChainAtomPositions[residueIndex][atomName] = self.chain[residueIndex].getAtom(atomName).getPosition()

    def restoreChain(self):
        for residueIndex in self.chain.residueRange():
            for atomName in self.chain[residueIndex].getAtomNames():
                self.chain[residueIndex].getAtom(atomName).setPosition(self.backupChainAtomPositions[residueIndex][atomName])

    def getClosestVertex(self, vertex, otherVertices, threshold, oneToOne):
        closestVertexIndex    = None
        closestVertex         = None
        closestVertexDistance = None
        for otherVertexIndex, otherVertex in enumerate(otherVertices):
            vertexDistance = (vertex - otherVertex).length()

            if (threshold == None or vertexDistance <= threshold) and (closestVertexDistance == None or vertexDistance < closestVertexDistance):
                closestVertexIndex    = otherVertexIndex
                closestVertex         = otherVertex
                closestVertexDistance = vertexDistance

        if oneToOne and closestVertexIndex != None:
            otherVertices.pop(closestVertexIndex)

        return closestVertex

    def calculateRMSD(self, sourceVertices, targetVertices, threshold, oneToOne):
        squaredDistanceSum = 0
        numSourceVertices  = len(sourceVertices)
        for sourceVertex in sourceVertices:
            targetVertex = self.getClosestVertex(sourceVertex, targetVertices, threshold, oneToOne)

            squaredDistanceSum += ((sourceVertex - targetVertex).length() ** 2)

        rmsd = sqrt(squaredDistanceSum / numSourceVertices)

        return rmsd

    def ntuples(self, lst, n):
        tuples = zip(*[lst[i:]+lst[:i] for i in range(n)])
        tuples.pop()

        return tuples

    def computeEuclideanDistance(p1, p2):
        # The Euclidean distance between two n dime points. Two input should have the same dim.
        distance = 0
        for i in range(0, len(p1)):
            distance += (p1[i] - p2[i])**2

        return distance**0.5

    def computeMidPoint(p1, p2):
        # The midpoint of p1 and p2.
        return [sum(x)*0.5 for x in zip(p1, p2)]
