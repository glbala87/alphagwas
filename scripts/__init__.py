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
