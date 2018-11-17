#!/usr/bin/env python
import sys
import time
import random
from functools import wraps


# Retry a function a certain number of times before raising an exception
def retry(ExceptionToCheck, tries=4, delay=4, backoff=1.5, logger=None):

    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):

            mtries, mdelay = tries, delay

            while mtries > 1:

                try:
                    return f(*args, **kwargs)

                except ExceptionToCheck:
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff

            return f(*args, **kwargs)

        return f_retry

    return deco_retry

# Print a status message when a function is running
class printstatus(object):

    def __init__(self, status):
        self.status = status

    def __call__(self, f, *args, **kwargs):

        def new_func(*args, **kwargs):
            sys.stdout.write(self.status + '... ')
            sys.stdout.flush()

            func = f(*args, **kwargs)

            sys.stdout.write('Done.\n')
            sys.stdout.flush()

            return func

        return new_func
