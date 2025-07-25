# -*- coding: utf-8 -*-

"""
/***************************************************************************
 Syn_SC
                                 A QGIS plugin
 Generate high-volume synthetic point data with user-defined continuity & smoothness
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2025-06-23
        copyright            : (C) 2025 by Anonymous Research Team
        email                : withheld@review.local
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Anonymous Research Team'
__date__ = '2025-06-23'
__copyright__ = '(C) 2025 by Anonymous Research Team'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from qgis.core import QgsProcessingProvider
#from .syn_sc_algorithm import Syn_SCAlgorithm
from .synsc_bruteforce import SynSCBrute
from .synsc_heuristic import SynSCHeuristic
from .synsc_iterative import SynSCIterative
from .synsc_scale_assistant import SynSCScaleAssistant
from .synsc_metric_evaluator import SynSCMetricEvaluator


class Syn_SCProvider(QgsProcessingProvider):

    def __init__(self):
        """
        Default constructor.
        """
        QgsProcessingProvider.__init__(self)

    def unload(self):
        """
        Unloads the provider. Any tear-down steps required by the provider
        should be implemented here.
        """
        pass

    def loadAlgorithms(self):
        """
        Loads all algorithms belonging to this provider.
        """
        #self.addAlgorithm(Syn_SCAlgorithm())
        self.addAlgorithm(SynSCBrute())
        self.addAlgorithm(SynSCHeuristic())
        self.addAlgorithm(SynSCIterative())
        self.addAlgorithm(SynSCScaleAssistant())
        self.addAlgorithm(SynSCMetricEvaluator())
        # add additional algorithms here
        # self.addAlgorithm(MyOtherAlgorithm())

    def id(self):
        """
        Returns the unique provider id, used for identifying the provider. This
        string should be a unique, short, character only string, eg "qgis" or
        "gdal". This string should not be localised.
        """
        return 'Syn-SC Toolkit'

    def name(self):
        """
        Returns the provider name, which is used to describe the provider
        within the GUI.

        This string should be short (e.g. "Lastools") and localised.
        """
        return self.tr('Syn-SC Toolkit')

    def icon(self):
        """
        Should return a QIcon which is used for your provider inside
        the Processing toolbox.
        """
        return QgsProcessingProvider.icon(self)

    def longName(self):
        """
        Returns the a longer version of the provider name, which can include
        extra details such as version numbers. E.g. "Lastools LIDAR tools
        (version 2.2.1)". This string should be localised. The default
        implementation returns the same string as name().
        """
        return self.name()
