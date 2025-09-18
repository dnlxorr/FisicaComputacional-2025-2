from abc import ABC, abstractmethod
import numpy as np

class Particle (ABC):
    def __init__(self, mass=0.0, charge=0.0, position=np.array([0.0, 0.0]), velocity = np.array([0.0, 0.0]), acceleration = np.array([0.0, 0.0])):
        self.mass=mass
        self.charge=charge
        self.position= position
        self.velocity = velocity
        self.acceleration = acceleration

        @abstractmethod
        def move(self):
            pass

        @abstractmethod
        def setPosition(self):
            pass

        @abstractmethod
        def getPosition(self):
            pass

        @abstractmethod
        def getVelocity(self):
            pass

        @abstractmethod
        def getAcceleration(self):
            pass


class Electron(Particle):

    def __init__(self):
        super().__init__(mass=9.11e-31,charge=-1.60e-19)

    def setPosition(self,pos):
        self.position = pos

    def getPosition(self):
        return self.position

    def setVelocity(self,vel):
        self.velocity = vel

    def getVelocity(self):
        return self.velocity

    def setAcceleration(self,acel):
         self.acceleration = acel

    def getAcceleration(self):
        return self.acceleration

    def move(self):
        self.position = self.position + np.array([1,1])

    def emisionFoton(self,energia):
        print(f"El electrón emite un fotón de {energia}")

    def calDeBroglieLongitudOnda(self):
        h = 6.62607015e-34
        v=np.linalg.norm(self.velocity) #Calcular el modulo

        if v==0:
            return mp.inf
        return h/(self.mass*v)


class Proton(Particle):

    def __init__(self):
        super().__init__(mass=1.6726e-27, charge=1.60e-19)

    def setPosition(self, pos):
        self.position = pos

    def getPosition(self):
        return self.position

    def setVelocity(self, vel):
        self.velocity = vel

    def getVelocity(self):
        return self.velocity

    def setAcceleration(self, acel):
        self.acceleration = acel

    def getAcceleration(self):
        return self.acceleration

    def move(self):
        self.position = self.position + np.array([1, 1])

    def calcularMomentoMagnetico(self):
        magneton_nuclear = 5.050783699e-27
        g_factor = 5.5856946893
        return g_factor*magneton_nuclear

class Neutron(Particle):

    def __init__(self):
        super().__init__(mass=1.6749e-27, charge=0)

    def setPosition(self, pos):
        self.position = pos

    def getPosition(self):
        return self.position

    def setVelocity(self, vel):
        self.velocity = vel

    def getVelocity(self):
        return self.velocity

    def setAcceleration(self, acel):
        self.acceleration = acel

    def getAcceleration(self):
        return self.acceleration

    def move(self,dt):
        self.velocity += self.acceleration * dt
        self.position += self.velocity * dt

    def calMomentoMagnetico(self):
        magneton_nuclear = 5.050783699e-27
        g_factor = -3.82608545
        return g_factor*magneton_nuclear
