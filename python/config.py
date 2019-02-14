"""This module supports instance configuration for so3g.  Eventually
this might be loaded from config file, perhaps ~/.so3g.  But for now it's
hard-coded to return sensible defaults.

"""
DEFAULTS = {
    'patch_g3frame': True,  # boolean.
    'use_astropy': 'try',   # string.  yes, no, try
    'use_pixell': 'try',    # string.  yes, no, try
}

def get_config():
    return DEFAULTS.copy()

