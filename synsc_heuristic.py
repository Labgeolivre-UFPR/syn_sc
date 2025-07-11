# -*- coding: utf-8 -*-
"""
QGIS Processing algorithm: Syn-SC Heuristic Generator
Implements the checker-board heuristic for minimum smoothness.

Author: <your name>, 2025
Licence: MIT
"""

from itertools import product
import random

# --- imports originais -------------------------------------------------
from qgis.core import QgsProcessingAlgorithm, QgsProcessingParameterVectorLayer
from .synsc_metric_evaluator import calculate_global_smoothness, neighbour_table
import random



from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsProcessingOutputString,
    QgsField,
    QgsFields,
    QgsFeature,
    QgsFeatureSink, 
    QgsGeometry,
    QgsPointXY,
    QgsSpatialIndex,
    QgsWkbTypes,
    QgsProcessingException
)
@staticmethod  
def tr(txt):
    return QCoreApplication.translate('SynSC', txt)
    return QCoreApplication.translate('SynSCHeuristic', message)

# ----------------------------------------------------------------------
def checkerboard_min(c_tgt: float, layer) -> float:
    """
    Returns the theoretical minimum smoothness (SM_min) for a grid that already has
    has the target continuity c_tgt. Purely mathematical operation;
    does not create a layer on disc.

    Parameters
    ----------
    c_tgt Percentage of populated cells (0-100).
    layer: hexagon grid layer already clipped to AOI.
    """
    feat_ids = [f.id() for f in layer.getFeatures()]
    n_total  = len(feat_ids)
    n_pop   = round(n_total * c_tgt / 100)
    neigh   = neighbour_table(layer)
    n_empty = n_total - n_pop

    # mapping ID ➜ index 0…n_total-1
    id2idx = {fid: idx for idx, fid in enumerate(feat_ids)}

    # ---------- Phase 1: checkerboard pattern in populated cells --------------------
    vals = [(100 if (i % 2 == 0) else 1) for i in range(n_pop)]
    vals.extend([0] * (n_total - n_pop))

    raw_neigh = neighbour_table(layer)            # keys = feature IDs
    neigh = {
        id2idx[fid]: [id2idx[nid] for nid in raw_neigh[fid] if nid in id2idx]
        for fid in raw_neigh
    }

    # ---------- Phase 2: high-contrast reinforcement --------------------
    # populated cells with ≥ 50 % empty neighbours → 100
    for i in range(n_pop):
        empty_neigh = sum(1 for j in neigh[i] if vals[j] == 0)
        if empty_neigh >= len(neigh[i]) / 2:
            vals[i] = 100

    # ---------- Phase 3: low contrast correction -------------------
    for i in range(n_pop):
        pop_neigh = [j for j in neigh[i] if vals[j] > 0]
        if pop_neigh and all(vals[j] == vals[i] for j in pop_neigh):
            # vira metade para o valor oposto
            for j in pop_neigh[::2]:
                vals[j] = 1 if vals[i] == 100 else 100

    # ---------- Global smoothness --------------------------------------
    sm_min = calculate_global_smoothness(vals, neigh)
    return sm_min
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def checkerboard_max(c_tgt: float, layer) -> float:
    """
    Upper bound of global smoothness (%) for a hex-grid that already
    fulfils the target continuity c_tgt.  Works in O(n).

    Strategy: give every populated cell the *opposite* extreme (1 or 100)
    of its axial-parity so that the largest possible proportion of
    populated–populated edges have |99| contrast.
    """
    # --- 1  IDs ↔ contiguous indices -----------------------------------
    feat_ids = [f.id() for f in layer.getFeatures()]
    id2idx   = {fid: i for i, fid in enumerate(feat_ids)}
    n_total  = len(feat_ids)
    
    # --- 2  neighbour list (indices) -----------------------------------
    raw_ngh = neighbour_table(layer)                         # IDs
    ngh = { id2idx[f]: [id2idx[n] for n in raw_ngh[f]]
            for f in raw_ngh }

    # --- 3  decide which cells stay empty ------------------------------
    n_pop   = round(n_total * c_tgt / 100)
    pop_ids = feat_ids[:n_pop]            # first n_pop = populated
    empty   = set(feat_ids[n_pop:])

    # --- 4  checker-board assignment -----------------------------------
    vals = [0] * n_total
    for fid in pop_ids:
        i = id2idx[fid]
        # axial parity: even = 1, odd = 100   (or swap – it’s symmetrical)
        vals[i] = 1  if (i % 2 == 0) else 100

    # --- 5  compute global SM ------------------------------------------
    sm = calculate_global_smoothness(vals, ngh)
    return sm
# ----------------------------------------------------------------------


class SynSCHeuristic(QgsProcessingAlgorithm):

    INPUT_GRID    = 'INPUT_GRID'
    C_TGT         = 'C_TGT'          
    P_MAX         = 'P_MAX'
    OUTPUT_POINTS = 'OUTPUT_POINTS'

    # ------------------------------------------------------------------ #
    # 1  GUI parameters
    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_GRID,
                tr('Hexagon grid'),
                [QgsProcessing.TypeVectorPolygon]))
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.C_TGT,
                tr('Continuity target (%)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=100, minValue=0, maxValue=100))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.P_MAX,
                tr('Maximum points per populated cell (p_max)'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=100, minValue=1, maxValue=100))

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_POINTS,
                tr('Synthetic point layer')))

    # ------------------------------------------------------------------ #
    # 2  Core processing routine
    def processAlgorithm(self, params, context, feedback):

        grid = self.parameterAsVectorLayer(params, self.INPUT_GRID, context)
        if grid is None:
            raise QgsProcessingException(tr('Input grid not found.'))
        c_tgt = self.parameterAsDouble(params, self.C_TGT,   context)
        p_max = self.parameterAsInt(params, self.P_MAX, context)
        fid_list = [f.id() for f in grid.getFeatures()]
        n_total  = len(fid_list)

        # ------- Prepare neighbour list using spatial index -------------
        idx = QgsSpatialIndex(grid.getFeatures())
        neighbours = {}
        empties = []          # cells already 0 (continuity initialised)
        for f in grid.getFeatures():
            geom = f.geometry()
            neigh_ids = [cid for cid in idx.intersects(geom.boundingBox())
                         if cid != f.id()
                         and geom.touches(grid.getFeature(cid).geometry())]
            neighbours[f.id()] = neigh_ids
        # ------------------------------------------------------
        # Build neighbour list and optionally detect empty cells
        idx = QgsSpatialIndex(grid.getFeatures())
        has_value_field = 'value' in [f.name() for f in grid.fields()]
        neighbours = {}
        empties = []

        for f in grid.getFeatures():
            geom = f.geometry()
            neigh_ids = [cid for cid in idx.intersects(geom.boundingBox())
                         if cid != f.id()
                         and geom.touches(grid.getFeature(cid).geometry())]
            neighbours[f.id()] = neigh_ids
            if has_value_field and f['value'] == 0:
                empties.append(f.id())
        # ------------------------------------------------------
        sm_min = checkerboard_min(c_tgt, grid)
        # ------------ Step 1: initial checkerboard ----------------------
        values = {}
        for f in grid.getFeatures():
            if f.id() in empties:
                values[f.id()] = 0
            else:
                # checkerboard: (row + col) parity
                centroid = f.geometry().centroid().asPoint()
                parity = (int(round(centroid.x()/10)) +
                          int(round(centroid.y()/10))) % 2
                values[f.id()] = 100 if parity == 0 else 1

        # ------------ Step 2: high-contrast reinforcement ---------------
        for fid, v in values.items():
            if v == 0:
                continue
            neighs = neighbours[fid]
            empty_count = sum(1 for n in neighs if values[n] == 0)
            if empty_count >= len(neighs) / 2:
                values[fid] = 100  # force to high value

        # ------------ Step 3: low-contrast correction -------------------
        for fid, v in values.items():
            if v == 0:
                continue
            neighs = [n for n in neighbours[fid] if values[n] > 0]
            if not neighs:
                continue
            same_count = sum(1 for n in neighs if values[n] == v)
            if same_count >= len(neighs) / 2:
                # flip half of those neighbours
                to_flip = neighs[: same_count // 2]
                for n in to_flip:
                    values[n] = 1 if values[n] == 100 else 100

        # ------------ Rescale to p_max (Eq. 3) --------------------------
        for fid, v in values.items():
            values[fid] = round(v * p_max / 100)

        # ------------ Create output sink --------------------------------
        out_fields = QgsFields()
        out_fields.append(QgsField('value', QVariant.Int))

        sink, dest_id = self.parameterAsSink(
            params, self.OUTPUT_POINTS, context,
            out_fields, QgsWkbTypes.Point, grid.sourceCrs())
        if sink is None:
            raise QgsProcessingException(tr('Failed to create output layer.'))

        # ------------ Scatter points (bounding-box rejection) -----------
        random.seed()
        for f in grid.getFeatures():
            count = values[f.id()]
            if count == 0:
                continue
            geom = f.geometry()
            if geom.isMultipart():
                ring = geom.asMultiPolygon()[0][0]      
            else:
                ring = geom.asPolygon()[0]              
            xs = [pt.x() for pt in ring]
            ys = [pt.y() for pt in ring]
            xmin, xmax = min(xs), max(xs)
            ymin, ymax = min(ys), max(ys)

            for _ in range(count):
                while True:
                    x = random.uniform(xmin, xmax)
                    y = random.uniform(ymin, ymax)
                    if geom.contains(QgsGeometry.fromPointXY(QgsPointXY(x, y))):
                        feat = QgsFeature(out_fields)
                        feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
                        feat['value'] = count
                        sink.addFeature(feat, QgsFeatureSink.FastInsert)
                        break

        return {self.OUTPUT_POINTS: dest_id}

    # ------------------------------------------------------------------ #
    # 3  Metadata
    def name(self):
        return 'synsc_heuristic'

    def displayName(self):
        return tr('Syn-SC Heuristic Generator')

    def group(self):
        return tr('Syn-SC')

    def groupId(self):
        return 'synsc'

    def createInstance(self):
        return SynSCHeuristic()

    def shortHelpString(self):
        return tr(
            'Generates a synthetic point layer using the checker-board '
            'heuristic that approximates the minimum global smoothness '
            'given the current continuity pattern.\n\n'
            'Recommended when n_total > 16 or target smoothness ≤ 40 %.')
