# -*- coding: utf-8 -*-
"""
QGIS Processing algorithm: Syn-SC Iterative Generator
Implements the ±(1–10) optimiser for target smoothness.

Author: <your name>, 2025
Licence: MIT
"""

import random
from math import fabs

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
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

def tr(txt):
    return QCoreApplication.translate('SynSC', txt)


class SynSCIterative(QgsProcessingAlgorithm):

    INPUT_GRID    = 'INPUT_GRID'
    C_TGT         = 'C_TGT'
    S_TGT         = 'S_TGT'
    P_MAX         = 'P_MAX'
    MAX_ITER      = 'MAX_ITER'
    OUTPUT_POINTS = 'OUTPUT_POINTS'

    # ------------------------------------------------------------------ #
    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_GRID,
                tr('Hexagon grid (must contain “value” field)'),
                [QgsProcessing.TypeVectorPolygon]))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.C_TGT,
                tr('Continuity target – populated cells (%)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=100, minValue=0, maxValue=100))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.S_TGT,
                tr('Smoothness target (%)'),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=25, minValue=0, maxValue=100))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.P_MAX,
                tr('Maximum points per populated cell (p_max)'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=100, minValue=1, maxValue=100))

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_ITER,
                tr('Maximum iterations'),
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=10000, minValue=1))

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_POINTS,
                tr('Synthetic point layer')))

    # ------------------------------------------------------------------ #
    def processAlgorithm(self, params, context, feedback):

        grid = self.parameterAsVectorLayer(params, self.INPUT_GRID, context)
        if grid is None:
            raise QgsProcessingException(tr('Input grid not found.'))

        c_tgt  = self.parameterAsDouble(params, self.C_TGT, context)
        s_tgt  = self.parameterAsDouble(params, self.S_TGT, context)
        p_max  = self.parameterAsInt(params, self.P_MAX, context)
        max_it = self.parameterAsInt(params, self.MAX_ITER, context)

        fid_list = [f.id() for f in grid.getFeatures()]
        n_total  = len(fid_list)
        n_empty_req = round(n_total * (100 - c_tgt) / 100)

        # ---------------- Build neighbour list --------------------------
        idx = QgsSpatialIndex(grid.getFeatures())
        neighbours = {}
        for f in grid.getFeatures():
            geom = f.geometry()
            neigh_ids = [cid for cid in idx.intersects(geom.boundingBox())
                         if cid != f.id()
                         and geom.touches(grid.getFeature(cid).geometry())]
            neighbours[f.id()] = neigh_ids

        # ---------------- Continuity initialisation --------------------
        # 1. stronger initial contrast
        values = {
            fid: (1 if i % 2 else 100)
            for i, fid in enumerate(fid_list)
        }
        if n_empty_req > 0:
            empty_ids = random.sample(fid_list, n_empty_req)
            for eid in empty_ids:
                values[eid] = 0
        # ---------------- Iterative optimiser --------------------------
        # ---------------------------------------------------------------
        #  Smoothness helper
        # ---------------------------------------------------------------
        def global_sm(vals):
            diffs = []
            for fid, v in vals.items():
                if v == 0:
                    continue
                neigh_vals = [vals[n] for n in neighbours[fid] if vals[n] > 0]
                diffs.extend(abs(v - nv) / 99 for nv in neigh_vals)
            return 100 * sum(diffs) / len(diffs) if diffs else 0

        # ---------------------------------------------------------------
        #  Iterative optimisation loop (accept-if-better)
        # ---------------------------------------------------------------
        sm = global_sm(values)
        feedback.pushInfo(
            tr(f'Initial SM = {sm:.2f} % (target {s_tgt:.2f})'))

        it = 0
        while fabs(sm - s_tgt) > 0.01 * s_tgt and it < max_it:

            # --- choose a random populated cell ------------------------
            populated = [fid for fid, v in values.items() if v > 0]
            if not populated:
                break
            fid = random.choice(populated)

            neigh_vals = [values[n] for n in neighbours[fid] if values[n] > 0]
            if not neigh_vals:
                it += 1
                continue
            neigh_mean = sum(neigh_vals) / len(neigh_vals)

            # --- adaptive step size (1 … 99) ---------------------------
            gap   = abs(sm - s_tgt) / s_tgt          # 0 … 1
            delta = max(1, int(gap * 99))            # larger when far away

            old_val = values[fid]

            # --- 1. propose new value (direction heuristic) ------------
            if sm < s_tgt:                           # need MORE contrast
                new_val = old_val + delta if old_val <= neigh_mean else old_val - delta
            else:                                    # need LESS contrast
                new_val = old_val - delta if old_val >= neigh_mean else old_val + delta
            new_val = max(1, min(100, new_val))      # clamp to 1–100

            # --- 2. evaluate global SM with the change -----------------
            values[fid] = new_val                    # tentative
            new_sm = global_sm(values)

            # --- 3. accept only if closer to the target ----------------
            if abs(new_sm - s_tgt) < abs(sm - s_tgt):
                sm = new_sm                          # improvement → keep
            else:
                values[fid] = old_val                # revert

            # --- bookkeeping & progress bar ---------------------------
            it += 1
            if it % 500 == 0:
                feedback.setProgress(int(100 * it / max_it))
                feedback.pushInfo(tr(f'Iter {it}: SM = {sm:.2f} %'))

        feedback.pushInfo(
            tr(f'Finished after {it} iterations: SM = {sm:.2f} %'))


        # ---------------- Rescale to p_max (Eq. 3) ----------------------
        for fid, v in values.items():
            values[fid] = round(v * p_max / 100)

        # ---------------- Create sink & scatter points -----------------
        out_fields = QgsFields()
        out_fields.append(QgsField('value', QVariant.Int))

        sink, dest_id = self.parameterAsSink(
            params, self.OUTPUT_POINTS, context,
            out_fields, QgsWkbTypes.Point, grid.sourceCrs())
        if sink is None:
            raise QgsProcessingException(tr('Failed to create output layer.'))

        random.seed()
        for f in grid.getFeatures():
            count = values[f.id()]
            if count == 0:
                continue
            geom = f.geometry()
            # ---------------------------------------------------------------
            #  Choose the exterior ring, no matter if the geometry is single
            #  or multipart.  We only need its bounding box for rejection
            #  sampling, so any first exterior ring is fine.
            # ---------------------------------------------------------------
            if geom.isMultipart():
                # asMultiPolygon() → [[[ring0], [hole1], …], …]
                rings = geom.asMultiPolygon()
                if not rings:
                    continue          # skip geometry if for some reason it's empty
                outer_ring = rings[0][0]
            else:
                # single polygon → [[ring0], [hole1], …]
                rings = geom.asPolygon()
                if not rings:
                    continue
                outer_ring = rings[0]

            xs = [pt.x() for pt in outer_ring]
            ys = [pt.y() for pt in outer_ring]
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
    def name(self):
        return 'synsc_iterative'

    def displayName(self):
        return tr('Syn-SC Iterative Generator')

    def group(self):
        return tr('Syn-SC')

    def groupId(self):
        return 'synsc'

    def createInstance(self):
        return SynSCIterative()

    def shortHelpString(self):
        return tr(
            'Generates a synthetic point layer by iteratively adjusting '
            'cell values (± 1–10) until global smoothness matches the '
            'target within 1 %. Continuity target is applied first.')
