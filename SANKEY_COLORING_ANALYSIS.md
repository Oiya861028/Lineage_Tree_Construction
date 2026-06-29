# Sankey Plot Coloring Analysis - Moscot/Moslin Implementation

## Overview
This document summarizes how sankey plot coloring is handled in the Moscot library (used by Moslin) based on analysis of the codebase.

---

## 1. Sankey Function Hierarchy

### Main Plotting Function: `mtp.sankey()`
**Location:** `src/moscot/plotting/_plotting.py` (lines 127-197)

```python
def sankey(
    obj: Union[AnnData, "TemporalProblem", "LineageProblem", "SpatioTemporalProblem"],
    key: str = _constants.SANKEY,
    captions: Optional[List[str]] = None,
    title: Optional[str] = None,
    colors: Optional[Dict[str, float]] = None,
    alpha: float = 1.0,
    interpolate_color: bool = False,
    cmap: Union[str, mpl.colors.Colormap] = "viridis",
    ax: Optional[mpl.axes.Axes] = None,
    return_fig: bool = False,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    save: Optional[Union[str, pathlib.Path]] = None,
    **kwargs: Any,
) -> Optional[mpl.figure.Figure]:
```

### Data Preparation Function: `lp.sankey()`
**Location:** `src/moscot/problems/time/_mixins.py` (TemporalMixin class, lines 176-283)

```python
def sankey(
    self,
    source: K,
    target: K,
    source_groups: Str_Dict_t,
    target_groups: Str_Dict_t,
    threshold: Optional[float] = None,
    normalize: bool = False,
    forward: bool = True,
    restrict_to_existing: bool = True,
    order_annotations: Optional[Sequence[str]] = None,
    key_added: Optional[str] = _constants.SANKEY,
    **kwargs: Any,
) -> Optional[list[pd.DataFrame]]:
```

### Internal Rendering Function: `_sankey()`
**Location:** `src/moscot/plotting/_utils.py` (lines 52-201)

```python
def _sankey(
    adata: AnnData,
    key: str,
    transition_matrices: List[pd.DataFrame],
    captions: Optional[List[str]] = None,
    colorDict: Optional[Union[Dict[Any, str], mpl.colors.ListedColormap]] = None,
    title: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    ax: Optional[mpl.axes.Axes] = None,
    cont_cmap: Union[str, mpl.colors.Colormap] = "viridis",
    fontsize: float = 12.0,
    horizontal_space: float = 1.5,
    force_update_colors: bool = False,
    alpha: float = 1.0,
    interpolate_color: bool = False,
    side_bar_width: float = 0.02,
    **kwargs: Any,
) -> mpl.figure.Figure:
```

---

## 2. How Colors Are Assigned to Source (Left) Columns

### Color Assignment Logic

The source column colors come from **categorical annotations** stored in the AnnData object. Here's the process:

#### Step 1: Extract Color Dictionary
**Location:** `src/moscot/plotting/_utils.py`, lines 72-83

```python
if colorDict is None:
    # Ensure palette is set for the categorical key
    set_palette(adata=adata, key=key, cont_cmap=cont_cmap, force_update_colors=force_update_colors)
    
    # Build color dictionary from existing colors in adata.uns
    colorDict = {
        cat: adata.uns[f"{key}_colors"][i] 
        for i, cat in enumerate(adata.obs[key].cat.categories)
    }
else:
    # Verify all categories have colors
    missing = sorted(label for label in adata.obs[key].cat.categories if label not in colorDict)
    if missing:
        raise ValueError(f"The following labels have missing colors: `{missing}`.")
```

#### Step 2: Color Mapping in Source Columns
**Location:** `src/moscot/plotting/_utils.py`, lines 326-342 (`_get_cmap_norm()` function)

```python
def _get_cmap_norm(...):
    if row_annotation != "cell":
        # Map each category to its color
        row_color_dict = {
            row_adata.obs[row_annotation].cat.categories[i]: col
            for i, col in enumerate(row_adata.uns[f"{row_annotation}_colors"])
        }
        # Create list of colors for transition matrix index (source labels) - REVERSED order
        row_colors = [row_color_dict[cat] for cat in transition_matrix.index][::-1]
    else:
        row_colors = [0, 0, 0]  # Default black if no annotation
```

### Key Points about Source Column Coloring:

1. **Color Source:** `adata.uns["{annotation_key}_colors"]` 
   - These colors are typically assigned by scanpy when creating categorical annotations
   - Colors are in the same order as categories in `adata.obs[key].cat.categories`

2. **Unique Color Assignment:** 
   - **ONE unique color per unique category/label**
   - Not per row, but per category type
   - Example: All cells labeled "TypeA" get the same color

3. **Order:** Source colors are **REVERSED** (`[::-1]`)
   - This is because the sankey diagram displays from top to bottom

4. **No Parameter to Change:** 
   - Colors are automatically derived from the categorical annotation's color palette
   - User can pass custom `colors` dict to `mtp.sankey()` to override

---

## 3. Related Color Handling Functions

### For Interpolation Between Source and Target
**Location:** `src/moscot/plotting/_utils.py`, lines 511-518

```python
def _color_transition(c1: str, c2: str, num: int, alpha: float) -> List[str]:
    """Create interpolated colors between two colors."""
    # When interpolate_color=True, creates a gradient from source to target colors
```

### Column Color Creation
**Location:** `src/moscot/plotting/_utils.py`, lines 521-538

```python
def _create_col_colors(adata: AnnData, obs_col: str, subset: Union[str, List[str]]) -> Optional[mpl.colors.Colormap]:
    """Create colormap for target columns based on source category."""
    # Extracts color for a specific category and creates a gradient
```

---

## 4. `mtp.sankey()` Function Options

### Primary Parameters for Coloring:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `colors` | `Dict[str, float]` or `None` | `None` | Custom color dict mapping labels to colors. If `None`, uses `adata.uns["{key}_colors"]` |
| `cmap` | `str` or `mpl.colors.Colormap` | `"viridis"` | Colormap for heat map if using interpolation |
| `interpolate_color` | `bool` | `False` | If `True`, interpolates colors from source to target |
| `alpha` | `float` | `1.0` | Transparency (0=transparent, 1=opaque) |

### Example Usage from Notebooks:

**From tutorial.ipynb (line 316):**
```python
lp.sankey(
    source=8,
    target=12,
    source_groups={"state": order[::-1]},
    target_groups={"state": order[::-1]},
)
mtp.sankey(lp, figsize=(6, 5))
```

**From moslin_lp_analysis_notebook.ipynb (line 816):**
```python
lp.sankey(
    source=7.5,
    target=8,
    source_groups={"cell_type": ordered_cell_type},
    target_groups={"cell_type": ordered_cell_type},
)
mtp.sankey(lp, figsize=(10, 7))
```

---

## 5. Data Flow Summary

```
AnnData object
    ├─ adata.obs[key] (categorical column)
    │   └─ .cat.categories (unique values)
    └─ adata.uns[f"{key}_colors"] (color list)
           └─ Indexed by category position

    ↓

lp.sankey() computes transition matrices:
    └─ Returns list[pd.DataFrame] with rows/cols as category names

    ↓

mtp.sankey() plots the diagram:
    └─ Extracts colors from adata.uns["{key}_colors"]
    └─ Maps each category to its corresponding color
    └─ Applies to both source (left) and target (right) columns
    └─ Handles interpolation if requested
```

---

## 6. Key Findings

### How Source Column Coloring Works:

1. **Automatic Coloring:**
   - Colors come from categorical metadata stored in `adata.uns`
   - Each unique category/label gets ONE unique color
   - Not configurable via parameters for per-cell coloring

2. **Control Points:**
   - `source_groups` parameter in `lp.sankey()`: determines which categories are shown
   - `order_annotations` parameter: controls order (top to bottom)
   - `colors` parameter in `mtp.sankey()`: can override the color mapping entirely

3. **No Per-Row Unique Colors:**
   - All rows with the same category label share the same color
   - This is by design for categorical visualization

4. **Color Source Hierarchy:**
   1. If `colors` dict passed to `mtp.sankey()` → use that
   2. Else if `adata.uns["{key}_colors"]` exists → use that
   3. Else → use default colors from palette

### Parameters for Customization:

- **To change source column colors:** Modify `adata.uns["{key}_colors"]` before calling `mtp.sankey()`, or pass custom `colors` dict
- **To change which categories appear:** Use `source_groups` parameter in `lp.sankey()`
- **To reorder rows:** Use `order_annotations` parameter in `lp.sankey()`
- **To interpolate colors:** Set `interpolate_color=True` in `mtp.sankey()`

---

## 7. Related Functions

- `set_palette()`: Sets up color palette for categorical data
- `set_plotting_vars()`: Stores sankey data in `adata.uns['moscot_results']['sankey']`
- `_get_cmap_norm()`: Creates colormaps and norms for matplotlib
- `_sankey()`: Core rendering function

---

## 8. References

**Moscot GitHub:** https://github.com/theislab/moscot
**Moscot Documentation:** https://moscot.readthedocs.io/

**Key Files:**
- `src/moscot/plotting/_plotting.py` - Main sankey function
- `src/moscot/plotting/_utils.py` - Color handling utilities
- `src/moscot/problems/time/_mixins.py` - TemporalMixin.sankey() method
