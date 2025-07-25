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

import os
import sys
import inspect

from qgis.core import QgsProcessingAlgorithm, QgsApplication
from .syn_sc_provider import Syn_SCProvider

cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]

if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)


class Syn_SCPlugin(object):

    def __init__(self):
        self.provider = None

    def initProcessing(self):
        """Init Processing provider for QGIS >= 3.8."""
        self.provider = Syn_SCProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()

    def unload(self):
        QgsApplication.processingRegistry().removeProvider(self.provider)
