# -*- coding: utf-8 -*-
"""
QGIS Processing algorithm: Syn-SC Brute-force Generator
Creates a synthetic point layer whose smoothness is minimised
by exhaustive enumeration (only feasible for <= 16 hexagons).

Author: <your name>, 2025
Licence: MITQgsProcessingOutputString,
"""

from qgis.PyQt.QtCore import (
    QCoreApplication, QVariant
)
from qgis.PyQt.QtGui import QColor
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingOutputString,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsVectorLayer,
    QgsField,
    QgsFields,           
    QgsFeature,
    QgsFeatureSink,      
    QgsGeometry,
    QgsPointXY,         
    QgsProject,
    QgsSpatialIndex,
    QgsWkbTypes,
    QgsProcessingException,
 
)
from itertools import product, combinations
import random
from .synsc_metric_evaluator import neighbour_table, calculate_global_smoothness

__all__ = ['brute_force_min', 'SynSCBrute'] 


def brute_force_min(c_tgt: float, layer) -> float:
    """
    Min SM by exhaustive search (n_total ≤ 16).
    Populated cells can take 1 or 100.
    """
    ids      = [f.id() for f in layer.getFeatures()]
    n_total  = len(ids)
    n_pop    = round(n_total * c_tgt / 100)

    # id ↔ idx mapping
    id2idx   = {fid: i for i, fid in enumerate(ids)}
    neigh_id = neighbour_table(layer)
    neigh = {
        id2idx[f]: [id2idx[n] for n in neigh_id[f]] for f in neigh_id
    }

    best = 100.0
    idx_set = range(n_total)
    # 1) escolha de quais células ficam povoadas
    for pop_ids in combinations(idx_set, n_pop):
        empty_ids = set(idx_set) - set(pop_ids)
        # 2) para essas povoadas, enumere combinações (1 ou 100)
        for pattern in product((1, 100), repeat=n_pop):
            vals = [0]*n_total
            for i, idx in enumerate(pop_ids):
                vals[idx] = pattern[i]
            sm = calculate_global_smoothness(vals, neigh)
            if sm == 0:        
                continue       
            best = min(best, sm)
    return best

class SynSCBrute(QgsProcessingAlgorithm):
    @staticmethod
    def tr(msg: str) -> str:
        return QCoreApplication.translate('SynSCBrute', msg)
    
    INPUT_GRID   = 'INPUT_GRID'
    C_TGT        = 'C_TGT' 
    SM_MIN_TXT   = 'SM_MIN_TXT'

    def initAlgorithm(self, config=None):

        self.addParameter(QgsProcessingParameterVectorLayer(
            self.INPUT_GRID,
            self.tr('Hexagon grid (polygon layer)'),
            [QgsProcessing.TypeVectorPolygon]))

        self.addParameter(QgsProcessingParameterNumber(
            self.C_TGT,
            self.tr('Continuity target (%)'),
            QgsProcessingParameterNumber.Double,
            defaultValue=100, minValue=0, maxValue=100))
        self.addOutput(QgsProcessingOutputString(
            self.SM_MIN_TXT,
            self.tr('Minimum achievable smoothness (%)')))

    # ---------------------------------------------------------------------

    def processAlgorithm(self, params, context, feedback):
        grid   = self.parameterAsVectorLayer(params, self.INPUT_GRID, context)
        c_tgt  = self.parameterAsDouble(params, self.C_TGT,    context)
        n_total = grid.featureCount()

        # ----- sanity checks ------------------------------------------------
        if n_total > 16:
            raise QgsProcessingException(
                tr('Brute-force enumeration limited to 16 cells '
                   '(current grid has {0}).').format(n_total))

        sm_min = brute_force_min(c_tgt, grid)
        return {self.SM_MIN_TXT: f'{sm_min:.2f}'}
        

    # ---------------------------------------------------------------------

    def _global_smoothness(self, val_dict, neigh_dict):
        """Compute global SM (%) for a given value assignment."""
        diffs = []
        for fid, v in val_dict.items():
            neigh_vals = [
                val_dict[nid] for nid in neigh_dict[fid]
                if val_dict[nid] > 0 and v > 0]
            if not neigh_vals:
                continue
            diffs.extend(abs(v - nv) / 99 for nv in neigh_vals)
        return 100 * (sum(diffs) / len(diffs)) if diffs else 0

    # ---------------------------------------------------------------------

    def name(self):
        return 'synsc_brute'

    def displayName(self):
        return self.tr('Syn-SC Brute Generator')

    def group(self):
        return self.tr('Syn-SC')

    def groupId(self):
        return 'synsc'

    def createInstance(self):
        return SynSCBrute()

    def shortHelpString(self):
        return self.tr(
            'Exhaustive search of all {0,1} assignments to minimise global '
            'smoothness. Limited to grids ≤ 16 cells for tractability.\n\n'
            'Outputs a point layer whose number of points per cell equals\n'
            'round(value × p_max / 100).')
