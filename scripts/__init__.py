# AlphaGWAS Pipeline Scripts
from . import (
    extract_variants,
    get_ld_proxies,
    liftover,
    alphagenome_predict,
    score_variants,
    utils
)

# Optional modules - import if available
__all__ = [
    'extract_variants',
    'get_ld_proxies',
    'liftover',
    'alphagenome_predict',
    'score_variants',
    'utils'
]

# Static visualization module
try:
    from . import visualize
    __all__.append('visualize')
except ImportError:
    pass

# Interactive visualization module (requires plotly)
try:
    from . import visualize_interactive
    __all__.append('visualize_interactive')
except ImportError:
    pass

# Annotation module
try:
    from . import annotate
    __all__.append('annotate')
except ImportError:
    pass
