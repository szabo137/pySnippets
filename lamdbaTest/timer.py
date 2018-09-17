"""
a simple timing decorator
"""
import time
def timer(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__)
            kw['log_time'][name] = (te - ts)
        else:
            print '%r  %1.2e ' % \
                  (method.__name__, (te - ts))
        return result
    return timed

if __name__=='__main__':
    from timeit import timeit

    def testFunc(a,b):
        time.sleep(1e-4)
        return a*b

    @timer
    def testFuncTimer(**kw):
        return testFunc(1,2)

    logT = {}
    print testFuncTimer(log_time = logT)

    print "log: %s"%(str(logT))

    print testFuncTimer()

    def timeWrapper():
        return testFunc(1,2)

    print timeit(timeWrapper,number=1000)
