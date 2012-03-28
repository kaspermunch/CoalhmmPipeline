
class BufferedList:
    "Buffered Iterator for iterating over Data data list"
    def __init__(self, lst):
        self.lst = lst
        for i, j in zip(self.lst, sorted(self.lst, reverse=True)):
            assert i.start == j.start

    def __iter__(self):
        return self
    def next(self):
        data = False
        if len(self.lst):
            data = self.lst.pop()
            return data
        else:
            raise StopIteration
    def putback(self, data):
        self.lst.append(data)

        for i, j in zip(self.lst, sorted(self.lst, reverse=True)):
            if i.start != j.start:
                print "### error"
                print i
                print j
                assert 0


class BufferedIterator:
    "Buffered List for iterating over Data data input stream"
    def __init__(self, iterator):
        self.buffer = []
        self.iterator = iterator
    def __iter__(self):
        return self
    def next(self):
        data = False
        if len(self.buffer):
            data = self.buffer.pop()
        else:
            data = self.iterator.next()
        if data:
            return data
        else:
            raise StopIteration
    def putback(self, data):
        self.buffer.append(data)
