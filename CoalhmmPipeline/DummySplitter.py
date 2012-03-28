#Do nothing implementation of a splitter
#just an example of interface
#constructed using
#DummySplitter()

#splits mafs if they have too long stretches of Ns
#TODO: write it...
class DummySplitter:
	def __init__(self):
		pass
	
	#return a list of mafs representing the old mafs that have split
	#or a list containing just the original maf if no splitting is required	
	def split(self, maf):
		return [maf]
