class IntervalQuery():
    def __init__(self, notFound=None):
        self.data = []
        self.dirty = True
        self.notFound = notFound
        
    def put(self, start, stop, payload):
        self.data.append((start, stop, payload))
        self.dirty = True
        
    def query(self, x):
        if self.dirty:
            self.data.sort()
            self.dirty = False
            
        return self._contains(x, 0, len(self.data))
        
    #binary search            
    def _contains(self, x, start, end):
        if start >= end:
            return self.notFound

        m = start + (end - start)/2
    
        (a,b,payload) = self.data[m]
    
        if a <= x and x <= b:
            return payload
        
        if x < a:
            return self._contains(x, start, m)
        else:
            return self._contains(x, m+1, end)
