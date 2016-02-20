class StatPoint() :
    """
    This object is used to store the result of get_jump method of sweep_set_open
    """
    def __init__(self, A = "None") :
        if(A == "None") :
            self.value = "None"
            self.field = "None"
            self.up = "None"
            self.trace = "None"
            self.sweep_nbr = "None"
        else :
            self.value = A[0]
            self.field = A[1]
            self.up = A[2]
            self.trace = A[3]
            self.sweep_nbr = A[4]

    def PassArray(self):
        return [self.value, self.field, self.up, self.trace, self.sweep_nbr ]

    def InBetween(self, seuil1, seuil2) :
        if self.value > seuil1 and self.value < seuil2 :
            return True
        else :
            return False

    def Inferior(self, seuil):
        if self.value < seuil :
            return True
        else :
            return False

    def Superior(self, seuil):
        if self.value > seuil :
            return True
        else :
            return False
