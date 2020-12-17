"""
This package simply provides a universal way to import spt3g, either from a
bundled subpackage or from somewhere on the filesystem.
"""

spt3g = None

if spt3g is None:
    try:
        from . import spt3g
    except:
        # Not bundled
        import spt3g
