import glob
from os.path import dirname,basename,isfile
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
#导入目录下的模块 print(modules,__all__)
