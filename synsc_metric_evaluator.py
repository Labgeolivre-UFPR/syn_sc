# -*- coding: utf-8 -*-
"""
Syn-SC Metric Evaluator
Computes continuity (CON) and smoothness (SM) for a point layer
at a user-defined analysis scale or explicit hexagon diameter.

Author: <anonymous for review>, 2025   Licence: MIT
"""

import random
from math import fabs

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSink,
    QgsProcessingOutputString,
    QgsSpatialIndex,
    QgsWkbTypes,
    QgsFeatureSink,
    QgsProcessingException
)
from qgis import processing

def tr(txt):
    return QCoreApplication.translate('SynSC', txt)

# --- helper functions (top-level) -----------------------------------
def calculate_global_smoothness(values, neighbours):
    """
    Returns global smoothness (%) for an array 'values'
    (0 = empty cell) and a neighbour-list dict.
    """
    diffs = []
    for idx, v in enumerate(values):
        if v == 0:
            continue
        neigh_ids = neighbours[idx]
        for j in neigh_ids:
            vj = values[j]
            if vj > 0:
                diffs.append(abs(v - vj) / 99)
    return 100 * sum(diffs) / len(diffs) if diffs else 0

def neighbour_table(layer):
    """
    Return a dict {feat_id: [neigh_id, …]} for the given hexagon layer.
    Two cells are neighbours if their geometries touch.
    """
    idx = QgsSpatialIndex(layer.getFeatures())
    neigh = {}
    for f in layer.getFeatures():
        geom = f.geometry()
        neigh_ids = [cid for cid in idx.intersects(geom.boundingBox())
                     if cid != f.id()
                     and geom.touches(layer.getFeature(cid).geometry())]
        neigh[f.id()] = neigh_ids
    return neigh


class SynSCMetricEvaluator(QgsProcessingAlgorithm):

    INPUT_PTS   = 'INPUT_PTS'
    SCALE_DENOM = 'SCALE_DENOM'
    HEX_DIAM    = 'HEX_DIAM'
    RETURN_GRID = 'RETURN_GRID'
    OUTPUT_GRID = 'OUTPUT_GRID'
    CONT_TXT    = 'CONT_TXT'
    SMOOTH_TXT  = 'SMOOTH_TXT'

    # --------------------------------------------------------------
    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
            'AOI',
            tr('Area of Interest (polygon layer) – optional'),
            [QgsProcessing.TypeVectorPolygon],
            optional=True))

        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_PTS,
                tr('Point layer to evaluate'),
                [QgsProcessing.TypeVectorPoint]))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.SCALE_DENOM,
                tr('Map scale denominator (1 : n)'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=50000, minValue=1))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.HEX_DIAM,
                tr('Override hexagon diameter (m) – optional'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0, minValue=0))

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.RETURN_GRID,
                tr('Return the filled grid as output layer?'),
                defaultValue=False))

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_GRID,
                tr('Grid with point counts (optional)'),
                defaultValue='TEMPORARY_OUTPUT'))
        
        self.addParameter(
            QgsProcessingParameterNumber(
                'MIN_PTS',
                tr('Minimum points that define a populated cell'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=1,  # keep compatible
                minValue=1
        ))
        # ---------- initAlgorithm -----------------------------------------
        self.addOutput(
            QgsProcessingOutputString(          # use OutputString
                self.CONT_TXT,
                tr('Continuity (%)'))
        )

        self.addOutput(
            QgsProcessingOutputString(          # idem aqui
                self.SMOOTH_TXT,
                tr('Smoothness (%)'))
        )

    # --------------------------------------------------------------
    def processAlgorithm(self, params, context, feedback):

        pt_layer = self.parameterAsVectorLayer(params, self.INPUT_PTS, context)
        if pt_layer is None:
            raise QgsProcessingException(tr('Point layer not found.'))

        scale_denom = self.parameterAsInt(params, self.SCALE_DENOM, context)
        d_hex_in    = self.parameterAsDouble(params, self.HEX_DIAM, context)
        return_grid = self.parameterAsBool(params, self.RETURN_GRID, context)

        d_hex = d_hex_in if d_hex_in > 0 else 0.005 * scale_denom
        feedback.pushInfo(tr(f'Using hexagon diameter = {d_hex:.2f} m'))

        # 1) Create hex grid covering the point layer extent
        grid_res = processing.run('qgis:creategrid', {
            'TYPE': 4,                            # hexagon polygons
            'HSPACING': d_hex,
            'VSPACING': d_hex,
            'EXTENT': pt_layer,
            'CRS': pt_layer.crs(),
            'OUTPUT': 'memory:grid_hex'
        }, context=context, feedback=feedback)
        grid = grid_res['OUTPUT']

        aoi_layer = self.parameterAsVectorLayer(params, 'AOI', context)

        if aoi_layer:
            sel_res = processing.run('native:extractbylocation', {
                'INPUT': grid,
                'PREDICATE': [0],          # 0 = intersects
                'INTERSECT': aoi_layer,
                'OUTPUT': 'memory:grid_sel'
            }, context=context, feedback=feedback)
            grid = sel_res['OUTPUT']
            feedback.pushInfo(tr(f'Grid after AOI filter: {grid.featureCount()} cells'))


        # 2) Count points per cell
        count_res = processing.run('native:countpointsinpolygon', {
            'POLYGONS': grid,
            'POINTS': pt_layer,
            'FIELD': 'value',
            'OUTPUT': 'memory:grid_val'
        }, context=context, feedback=feedback)
        grid = count_res['OUTPUT']          # now has 'value' field
        min_pts = self.parameterAsInt(params, 'MIN_PTS', context)
        n_total = grid.featureCount()
        n_pop = sum(1 for f in grid.getFeatures() if f['value'] >= min_pts)
        cont    = 100 * n_pop / n_total if n_total else 0

        # 3) Build neighbour list
        idx = QgsSpatialIndex(grid.getFeatures())
        neighbours = {}
        neighbours = neighbour_table(grid)

        # 4) Smoothness
        def global_sm():
            diffs = []
            for f in grid.getFeatures():
                v = f['value']
                if v == 0:
                    continue
                neigh = [grid.getFeature(n) for n in neighbours[f.id()]
                         if grid.getFeature(n)['value'] > 0]
                diffs.extend(abs(v - n['value']) / 99 for n in neigh)
            return 100 * sum(diffs) / len(diffs) if diffs else 0

        sm = global_sm()

        feedback.pushInfo(tr(f'Continuity = {cont:.2f} %'))
        feedback.pushInfo(tr(f'Smoothness = {sm:.2f} %'))

        # 5) Optionally output the grid
        grid_id = None
        if return_grid:
            sink, grid_id = self.parameterAsSink(
                params, self.OUTPUT_GRID, context,
                grid.fields(), QgsWkbTypes.Polygon, grid.crs())
            if sink is None:
                raise QgsProcessingException(tr('Failed to create grid sink.'))
            for feat in grid.getFeatures():
                sink.addFeature(feat, QgsFeatureSink.FastInsert)

        return {
            self.OUTPUT_GRID: grid_id,
            self.CONT_TXT:   f'{cont:.2f}',
            self.SMOOTH_TXT: f'{sm:.2f}'
        }

    # ---------------- metadata -----------------
    def name(self):        return 'synsc_metric_evaluator'
    def displayName(self): return tr('Syn-SC Metric Evaluator')
    def group(self):       return tr('Syn-SC')
    def groupId(self):     return 'synsc'
    def createInstance(self): return SynSCMetricEvaluator()
    def shortHelpString(self):
        return tr(
            'Calculates continuity (percentage of populated cells) and '
            'smoothness (mean neighbour contrast) for a point layer. '
            'Hexagon size is either 0.005 × map scale denominator or an '
            'explicit diameter in metres.')
