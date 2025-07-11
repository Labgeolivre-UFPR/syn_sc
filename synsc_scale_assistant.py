# -*- coding: utf-8 -*-
"""
Syn-SC Grid Advisor  •  QGIS Processing algorithm
Calculates scale-dependent cell size and theoretical minimum smoothness.

Author: <your name>, 2025  •  Licence: MIT
"""

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterString,
    QgsProcessingOutputString,
    QgsFeatureSink,
    QgsWkbTypes,
    QgsProcessingException
)
from qgis import processing
from .synsc_heuristic   import checkerboard_min, checkerboard_max
from .synsc_bruteforce  import brute_force_min
from qgis.core          import QgsCoordinateTransformContext
import math

def tr(txt):
    return QCoreApplication.translate('SynSC', txt)


class SynSCScaleAssistant(QgsProcessingAlgorithm):
    @staticmethod
    def tr(msg: str) -> str:
        return QCoreApplication.translate('SynSCScaleAssistant', msg)
    AOI          = 'AOI'
    SCALE_DENOM  = 'SCALE_DENOM'
    DPI          = 'DPI'
    C_TGT          = 'C_TGT'
    OUTPUT_GRID  = 'OUTPUT_GRID'
    D_HEX_TXT    = 'D_HEX_TXT'
    SOLVER_TXT   = 'SOLVER_TXT'
    SM_FLOOR_TXT  = 'SM_FLOOR_TXT'
    SM_CEIL_TXT = 'SM_CEIL_TXT'

    # --------------------------------------------------------------
    def initAlgorithm(self, config=None):

        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.AOI,
                tr('Area of Interest (polygon layer)'),
                [QgsProcessing.TypeVectorPolygon]))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.SCALE_DENOM,
                tr('Map scale denominator (1 : n)'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=50000, minValue=1))
        
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_GRID,
                tr('Hexagonal grid (output)')))
        
        self.addParameter(
            QgsProcessingParameterNumber(
                'C_TGT',
                tr('Continuity target (%)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=60, minValue=0, maxValue=100))

        self.addOutput(
            QgsProcessingOutputString(
                self.D_HEX_TXT,
                self.tr('Recommended hexagon diameter (m)')))
        
        self.addOutput(
            QgsProcessingOutputString(
                self.SOLVER_TXT,
                self.tr('Floor with 1/100 values (%)')))
        
        self.addOutput(
            QgsProcessingOutputString(
                self.SM_FLOOR_TXT,
                self.tr('Minimum achievable smoothness (%)')))
        
        self.addOutput(QgsProcessingOutputString(
                self.SM_CEIL_TXT, self.tr('Ceiling with 1/100 values (%)')))


    # --------------------------------------------------------------
    def processAlgorithm(self, params, context, feedback):

        aoi_layer = self.parameterAsVectorLayer(params, self.AOI, context)
        if aoi_layer is None:
            raise QgsProcessingException(tr('AOI layer not found.'))

        scale_denom = self.parameterAsInt(params, self.SCALE_DENOM, context)
        dpi       = self.parameterAsInt   (params, self.DPI, context)
        c_tgt     = self.parameterAsDouble(params, self.C_TGT, context)

        # 1) recommended hexagon diameter (≈ 5 mm on the map)
        mm_screen = 5.0
        d_hex = mm_screen * scale_denom / 1000.0           # metres on the ground
        feedback.pushInfo(tr(f'Recommended d_hex = {d_hex:.2f} m'))

        # 2) create the grid using qgis:creategrid
        grid_res = processing.run('qgis:creategrid', {
            'TYPE': 4,                       # hexagon (polygon)
            'HSPACING': d_hex,
            'VSPACING': d_hex,
            'EXTENT': aoi_layer,
            'HOVERLAY': 0, 'VOVERLAY': 0,
            'CRS': aoi_layer.crs(),
            'OUTPUT': 'memory:grid_hex'
        }, context=context, feedback=feedback)

        grid = grid_res['OUTPUT']
        # 2a) clip grid to AOI  ← NEW
        select_res = processing.run(
            'native:extractbylocation', {
                'INPUT': grid,              # grade completa
                'PREDICATE': [0],           # 0 = intersects
                'INTERSECT': aoi_layer,
                'OUTPUT': 'memory:grid_sel'
            },
            context=context, feedback=feedback
        )
        grid = select_res['OUTPUT']
        n_total = grid.featureCount()
        feedback.pushInfo(tr(f'Grid after selection: {n_total} hexagons.'))
        sm_init = checkerboard_min(c_tgt, grid)
        
        # 3) decide which Syn-SC solver to suggest
        if n_total <= 16:
            solver = 'brute'
            sm_floor = brute_force_min(c_tgt, grid)
        elif n_total > 1000:
            solver = 'heuristic'
        else:
            solver = 'iterative'
       
        sm_floor = checkerboard_min(c_tgt, grid)
        sm_ceiling = checkerboard_max(c_tgt, grid)      # NEW
        feedback.pushInfo( self.tr(f'Smoothness range with 1/100 levels: '
            f'{sm_floor:.1f} % – {sm_ceiling:.1f} %'))
        feedback.pushInfo(tr(f'Suggested solver: {solver}'))
        feedback.pushInfo(tr('Iterative solver may produce smoother grids '
                    'because it is allowed to use any value 1 … 100.'))
        # 4) write grid to output sink
        sink, dest_id = self.parameterAsSink(
            params, self.OUTPUT_GRID, context,
            grid.fields(), QgsWkbTypes.Polygon, grid.crs())
        if sink is None:
            raise QgsProcessingException(tr('Failed to create output grid.'))

        for feat in grid.getFeatures():
            sink.addFeature(feat, QgsFeatureSink.FastInsert)

        # 5) return variables for use in larger models
        return {
            self.OUTPUT_GRID: dest_id,
            self.D_HEX_TXT:  f'{d_hex:.2f}',
            self.SOLVER_TXT: solver,
            self.SM_FLOOR_TXT: f'{sm_floor:.2f}',
            self.SM_CEIL_TXT:  f'{sm_ceiling:.2f}'
        }

    # ----------------- metadata ----------------------------------
    def name(self): return 'synsc_scale_assistant'
    def displayName(self): return tr('Syn-SC Scale Assistant')
    def group(self): return tr('Syn-SC')
    def groupId(self): return 'synsc'
    def createInstance(self): return SynSCScaleAssistant()
    def shortHelpString(self):
        return tr(
            'Suggests a hexagon diameter for a given map scale, generates the '
            'grid, counts cells, and recommends the Syn-SC generator '
            '(brute, heuristic or iterative) based on grid size and target '
            'smoothness.'
        )
