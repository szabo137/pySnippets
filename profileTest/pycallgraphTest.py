"""
test of PyCallGraph
"""

from pycallgraph import PyCallGraph, Config, GlobbingFilter
from pycallgraph.output import GraphvizOutput
config = Config()
config.trace_filter = GlobbingFilter(exclude=[
    'numpy.*',
    'scipy.*',
    'ctypes.*',
    'matplotlib.*',
    'unittest.*',
    'math.*',
    'numbers.*'
])




with PyCallGraph(output=GraphvizOutput(),config=config):
   import testMod
   testMod.runAll()
