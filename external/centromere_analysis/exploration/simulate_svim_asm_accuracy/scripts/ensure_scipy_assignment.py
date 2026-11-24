"""Recover transient missing-file reads without changing the installed environment."""
import importlib
import sys
import time


def ensure_assignment_available(attempts=5):
    for attempt in range(attempts):
        try:
            module=importlib.import_module('scipy.optimize')
            r,c=module.linear_sum_assignment([[0,2],[2,0]])
            if list(r)!=[0,1] or list(c)!=[0,1]:
                raise RuntimeError('Unexpected assignment solver preflight result')
            return
        except FileNotFoundError:
            if attempt+1==attempts:raise
            for name in list(sys.modules):
                if name=='scipy' or name.startswith('scipy.'):
                    sys.modules.pop(name,None)
            time.sleep(2)
