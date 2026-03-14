# AlphaGWAS Pipeline Scripts
from . import (
    extract_variants,
    get_ld_proxies,
    liftover,
    alphagenome_predict,
    score_variants,
    utils
)

# Core modules always available
__all__ = [
    'extract_variants',
    'get_ld_proxies',
    'liftover',
    'alphagenome_predict',
    'score_variants',
    'utils'
]

# Optional modules - import if dependencies available

# Static visualization (requires matplotlib)
try:
    from . import visualize
    __all__.append('visualize')
except ImportError:
    pass

# Interactive visualization (requires plotly)
try:
    from . import visualize_interactive
    __all__.append('visualize_interactive')
except ImportError:
    pass

# LocusZoom plots
try:
    from . import locuszoom
    __all__.append('locuszoom')
except ImportError:
    pass

# Variant annotation (requires requests)
try:
    from . import annotate
    __all__.append('annotate')
except ImportError:
    pass

# Fine-mapping integration
try:
    from . import finemapping
    __all__.append('finemapping')
except ImportError:
    pass

# Pathway enrichment
try:
    from . import enrichment
    __all__.append('enrichment')
except ImportError:
    pass

# Colocalization analysis
try:
    from . import colocalization
    __all__.append('colocalization')
except ImportError:
    pass

# PDF report generation (requires reportlab)
try:
    from . import report
    __all__.append('report')
except ImportError:
    pass

# Multi-phenotype comparison
try:
    from . import multiphenotype
    __all__.append('multiphenotype')
except ImportError:
    pass

# REST API (requires fastapi)
try:
    from . import api
    __all__.append('api')
except ImportError:
    pass
