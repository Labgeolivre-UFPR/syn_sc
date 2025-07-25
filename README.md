<h1>Syn-SC &ndash; Continuity&ndash;Smoothness-controlled Synthetic Point Data (QGIS 3)</h1>

<p>
<strong>Syn-SC</strong> is a QGIS Processing provider that generates large synthetic point
clouds while letting you specify <em>continuity</em>
(proportion of empty hexagons) and <em>smoothness</em>
(mean absolute contrast between neighbouring cells) independently.
A third parameter &ndash; the point-per-cell cap &ndash; disguises the synthetic
origin of the data so the output looks natural enough for perceptual
user studies, AI training and teaching.
</p>

<h2>Key features</h2>
<ul>
  <li>Scale-adaptive hexagon grid (5&nbsp;mm on map or 5&nbsp;px on screen)</li>
  <li>Two-value smoothness window tells you which targets will converge</li>
  <li>Three solvers<br>
    &bull;&nbsp;<code>Brute</code> &nbsp;exact floor for grids ≤ 16<br>
    &bull;&nbsp;<code>Iterative</code> converges rapidly inside the window<br>
    &bull;&nbsp;<code>Heuristic</code> reaches the floor on production-sized grids
  </li>
  <li>Metric evaluator checks continuity &amp; smoothness post hoc</li>
  <li>Outputs Shapefile or GeoPackage; CRS preserved</li>
</ul>

<h2>Installation</h2>
<ol>
  <li>Download the latest ZIP from the GitHub
      <a href="releases">releases</a>.</li>
  <li>In QGIS &mdash; <em>Plugins &rsaquo; Manage and Install &rsaquo; Install from ZIP</em>.</li>
  <li>Restart QGIS, open the <em>Processing&nbsp;Toolbox</em> &nbsp;&rarr;&nbsp; <code>Syn-SC</code>.</li>
</ol>

<h2>Quick start</h2>
<ol>
<li>Run <strong>Scale Assistant</strong> on an AOI polygon; note the suggested solver and the two-value smoothness window.</li>
<li>Execute the recommended solver (Brute / Iterative / Heuristic).</li>
<li>(Optional) run <strong>Metric Evaluator</strong> to confirm that continuity and smoothness match your targets.</li>
</ol>

<h2>Parameters at a glance</h2>
<table>
<thead><tr><th>Symbol</th><th>Name&nbsp;/ Function</th><th>Typical range</th></tr></thead>
<tbody>
<tr><td><code>c<sub>tgt</sub></code></td><td>Target continuity (%)</td><td>10–90</td></tr>
<tr><td><code>s<sub>tgt</sub></code></td><td>Target smoothness (%)</td><td>0–100*</td></tr>
<tr><td><code>d<sub>hex</sub></code></td><td>Hexagon diameter (m)</td><td>50&nbsp;–&nbsp;2000</td></tr>
<tr><td><code>p<sub>max</sub></code></td><td>Point cap per cell</td><td>10&nbsp;–&nbsp;100</td></tr>
</tbody>
</table>
<p><small>*targets below the binary floor or above the ceiling are rejected unless the Iterative solver can satisfy them.</small></p>

<h2>Algorithmic details</h2>
<ol>
  <li><strong>Scale Assistant</strong> &ndash; builds grid, computes binary floor (SM<sub>floor</sub>) and ceiling (SM<sub>ceil</sub>), suggests solver.</li>
  <li><strong>Brute</strong> / <strong>Heuristic</strong> / <strong>Iterative</strong> – assign values 0,1&hellip;100 to meet <code>c<sub>tgt</sub></code> and <code>s<sub>tgt</sub></code>.</li>
  <li>Point realisation – rescale <code>value</code>&nbsp;&rarr; point count and scatter points uniformly.</li>
</ol>

<h2>Citing</h2>
<p>
If Syn-SC is helpful in your research, please cite:<br>
<em>Your&nbsp;Name (2025). Syn-SC: generating continuity- and smoothness-controlled synthetic point data.</em>
</p>

<h2>Licence</h2>
<p>MIT Licence – do anything, give credit, no warranty.</p>
