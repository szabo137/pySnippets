"""
testmod for cProfile
"""
import cProfile
import testMod
cProfile.run('testMod.runAll()',filename="myfile.profile")
