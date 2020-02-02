#!/usr/bin/sh
export PYTHONPATH=${PYTHONPATH}:/snap/inkscape/current/share/inkscape/extensions
/usr/bin/env python3 svx_output.py "$@"
