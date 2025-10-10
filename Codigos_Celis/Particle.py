from abc import ABC, abstractmethod
import math

# ==================================================
# Clase base abstracta
# ==================================================
class Particle(ABC):
    def __init__(self, mass=1.0, charge=0.0):
        self.mass = mass
        self.charge = charge
        self.pos = [0.0, 0.0]
        self.vel = [0.0, 0.0]
        self.acc = [0.0, 0.0]
        self.trayectoria_x = []
        self.trayectoria_y = []

    @abstractmethod
    def move(self, dt):
        pass

    @abstractmethod
    def applyBoundary(self, N):
        pass

    def setPosition(self, pos):
        self.pos = [pos[0], pos[1]]

    def getPosition(self):
        return [self.pos[0], self.pos[1]]

    def setVelocity(self, vel):
        self.vel = [vel[0], vel[1]]

    def getVelocity(self):
        return [self.vel[0], self.vel[1]]

    def setAcceleration(self, acel):
        self.acc = [acel[0], acel[1]]

    def getAcceleration(self):
        return [self.acc[0], self.acc[1]]

# ==================================================
# Clases concretas
# ==================================================

class IonTitanio(Particle):
    def __init__(self, x=0.0, y=0.0, vx=0.0, vy=0.0):
        super().__init__(mass=7.95e-26, charge=-1.602e-19)
        self.pos = [x, y]
        self.vel = [vx, vy]
        self.acc = [0.0, 0.0]
        self.trayectoria_x = [x]
        self.trayectoria_y = [y]

    def move(self, dt):
        # Movimiento tipo Leapfrog (simplificado)
        self.vel[0] += self.acc[0] * dt
        self.vel[1] += self.acc[1] * dt
        self.pos[0] += self.vel[0] * dt
        self.pos[1] += self.vel[1] * dt

        # Guardar trayectoria
        self.trayectoria_x.append(self.pos[0])
        self.trayectoria_y.append(self.pos[1])

    def applyBoundary(self, N):
        # Rebote en bordes (en lugar de desaparecer)
        if self.pos[0] > N:
            self.pos[0] = N
            self.vel[0] *= -1
        elif self.pos[0] < 0:
            self.pos[0] = 0
            self.vel[0] *= -1

        if self.pos[1] > N:
            self.pos[1] = N
            self.vel[1] *= -1
        elif self.pos[1] < 0:
            self.pos[1] = 0
            self.vel[1] *= -1


class Electron(Particle):
    def __init__(self, x=0.0, y=0.0, vx=0.0, vy=0.0):
        super().__init__(mass=9.11e-31, charge=-1.602e-19)
        self.pos = [x, y]
        self.vel = [vx, vy]
        self.acc = [0.0, 0.0]
        self.trayectoria_x = [x]
        self.trayectoria_y = [y]

    def move(self, dt):
        self.vel[0] += self.acc[0] * dt
        self.vel[1] += self.acc[1] * dt
        self.pos[0] += self.vel[0] * dt
        self.pos[1] += self.vel[1] * dt
        self.trayectoria_x.append(self.pos[0])
        self.trayectoria_y.append(self.pos[1])

    def applyBoundary(self, N):
        if self.pos[0] > N or self.pos[0] < 0:
            self.vel[0] *= -1
        if self.pos[1] > N or self.pos[1] < 0:
            self.vel[1] *= -1


class Proton(Particle):
    def __init__(self, x=0.0, y=0.0, vx=0.0, vy=0.0):
        super().__init__(mass=1.6726e-27, charge=1.602e-19)
        self.pos = [x, y]
        self.vel = [vx, vy]
        self.acc = [0.0, 0.0]
        self.trayectoria_x = [x]
        self.trayectoria_y = [y]

    def move(self, dt):
        self.vel[0] += self.acc[0] * dt
        self.vel[1] += self.acc[1] * dt
        self.pos[0] += self.vel[0] * dt
        self.pos[1] += self.vel[1] * dt
        self.trayectoria_x.append(self.pos[0])
        self.trayectoria_y.append(self.pos[1])

    def applyBoundary(self, N):
        if self.pos[0] > N or self.pos[0] < 0:
            self.vel[0] *= -1
        if self.pos[1] > N or self.pos[1] < 0:
            self.vel[1] *= -1
