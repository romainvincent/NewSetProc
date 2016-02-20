from numpy import linspace, size, linalg, real, deg2rad, cos, sin
from matplotlib.pyplot import plot
from setproc.mag.functions.mPauli import m_pauli



class SysMag():
    """
    This class allows to define easily a magnetic system with
    """
    def __init__(self, J, D, E, g):
        temp = m_pauli(J)
        self.Sx = temp[0]
        self.Sy = temp[1]
        self.Sz = temp[2]
        self.Sp = temp[3]
        self.Sm = temp[4]
        self.D = D
        self.E = E
        self.g = g

    def Hs(self, B, muB = 1,J = 0, theta = 0):
        TH = deg2rad(theta)
        Jx = J * cos(TH)
        Jz = J * sin(TH)
        H = - self.D * self.Sz**2
        H = H + Jz * self.Sz + Jx * self.Sx #introduce a coupling with a spin either up or down
        H = H + self.E * (self.Sp**2 + self.Sm**2)
        H = H + muB * (self.g[0] * self.Sx * B[0] + self.g[1] * self.Sy * B[1] + self.g[2] * self.Sz * B[2])
        return H

    def Zeeman(self, Bmin, Bmax, nbr, Bx = 0, By =0):
        self.zeem = dict()
        self.zeem["Bmin"] = Bmin
        self.zeem["Bmax"] = Bmax
        self.result = []
        B = linspace(Bmin, Bmax, nbr)
        for i in range(size(B)) :
            H =  self.Hs( [Bx , By, B[i]], 1)
            E = linalg.eigvals(H)
            self.result.append(real(E))
        return self.result

    def ZeemanPoint(self, Bx, By, Bz, muB = 1, J = 0, theta = 0) :
        H = self.Hs([Bx,By,Bz], muB, J, theta)
        E = linalg.eigvals(H)
        return real(E)



    def PlotZeeman(self, marker = 'b.'):
        R = dict()
        si = size(self.result[0])
        nbr = len(self.result)
        B = linspace(self.zeem["Bmin"], self.zeem["Bmax"], nbr)
        for j in range(si) :
            R[j] =[]

        for x in self.result :
            for j in range(si) :
                R[j].append(x[j])

        for j in range(si):
            plot(B, R[j], marker)

        return True

    def TbPc2()
        H_LF = -410/99*O_20 - 230*2/16335*O_40 - 30/891891*O_60