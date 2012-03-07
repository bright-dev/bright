class IOPair:
    coordinate = [0,0]
    left = coordinate[0]
    right = coordinate[1]

    def set_coordinate(self,x,y):
        self.coordinate = [x,y]
	self.left = self.coordinate[0]
	self.right = self.coordinate[1]

"""
pair = io_pair()
pair.set_coordinate(1,2)
print pair.coordinate, pair.left, pair.right
"""
print "Hello world"
