#Cell Class

class Cell(object):
    def __init__(self, typ, param):
        self.typ = 0
        self.fitness = param.GetFitness(self.typ)
        self.mutation = param.GetMutation(self.typ)
    def divide()
        
        
