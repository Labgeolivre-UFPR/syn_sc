**Purpose**  
Tessellates the Area of Interest with size-adaptive hexagons, computes the two-value smoothness window, and recommends the solver that will reach a user-supplied smoothness target in the shortest time.

**Parameters**  
| Name | Type | Default | Description |
|------|------|---------|-------------|
| *AOI layer* | Polygon | — | Spatial extent for synthesis. |
| *Map-scale denominator* | Integer | 50000 | Scale at which the output will be viewed. 5 mm on the map equals one cell diameter. |
| *Target continuity (%)* `c_tgt` | Double | 60 | Share of populated cells. |
| *Target smoothness (%)* `s_tgt` | Double | 50 | Desired global smoothness. |
| *DPI (optional)* | Integer | 96 | Use 5 px on-screen instead of 5 mm on paper. |

**Outputs**  
* Hexagon grid (polygon layer) with an empty **value** field  
* Text – recommended d_hex (m)  
* Text – two-value smoothness window (“floor – ceiling, %”)  
* Text – suggested solver (brute / iterative / heuristic)

**Algorithm**  
1. Derive cell diameter  
   \(d_{\text{hex}} = 5\text{ mm} \times \text{scale} / 1000\) (or 5 px / DPI).  
2. Tessellate AOI, build R-tree neighbour index.  
3. Compute SM_floor_2-level and SM_ceil_2-level under binary values {1, 100}.  
4. Recommend solver: Brute if *n*≤16; Iterative if 17–500 & target in window; otherwise Heuristic.

**Citations**  
MacEachren & DiBiase (1991) for continuity/smoothness concepts.  
*Your Name (2025) – Syn-SC plug-in*.
