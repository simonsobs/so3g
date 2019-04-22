import time

class Timer(): 
    def __init__(self, block_name=None): 
        self.t0 = time.time()
          
    def __enter__(self): 
        return self
      
    def report(self):
        return time.time() - self.t0

    def __exit__(self, exc_type, exc_value, exc_traceback): 
        print('Timer block exits after %.6f seconds' % (time.time() - self.t0))
