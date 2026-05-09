"""Entry point for ``python -m pbg_composite_nfsim_caspule``.

Defers to :func:`pbg_composite_nfsim_caspule.run.main`, which parses
the CLI args and runs the requested PBG document.
"""

import sys

from pbg_composite_nfsim_caspule.run import main

if __name__ == '__main__':
    sys.exit(main())
