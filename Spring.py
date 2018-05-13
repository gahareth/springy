import math
import matplotlib.pyplot as plt
import sys
from Core import AccelerationDueToGravity
from Core import RK4SecondOrder
from Core import Euler 
from Core import PhysicalQuantity
from Load import Load

class Spring:
    def __init__(self, InitialDisplacement, SpringConstant, Load, DampingCoefficient):
        self.Displacement = PhysicalQuantity(InitialDisplacement, 0)
        self.SpringConstant = SpringConstant
        self.Load = Load 
        self.DampingCoefficient = DampingCoefficient

    def Update(self, TimeStep):
        RK4SecondOrder(self.Displacement, TimeStep, lambda Time, Velocity, Displacement: self.Acceleration(Displacement, Velocity)) 

    def Acceleration(self, Displacement, Velocity):
        return (self.Load.Weight - self.DampingCoefficient * Velocity - self.SpringConstant * Displacement) / self.Load.Mass

class RungeKuttaSpring(Spring):
    def __init__(self, InitialDisplacement, SpringConstant, Load, DampingCoefficient):
        Spring.__init__(self, InitialDisplacement, SpringConstant, Load, DampingCoefficient)

    def Update(self, TimeStep):
        RK4SecondOrder(self.Displacement, TimeStep, lambda Time, Velocity, Displacement: self.Acceleration(Displacement, Velocity)) 

class EulerSpring(Spring):
    def __init__(self, InitialDisplacement, SpringConstant, Load, DampingCoefficient):
        Spring.__init__(self, InitialDisplacement, SpringConstant, Load, DampingCoefficient)

    def Update(self, TimeStep):
        Euler(self.Displacement, TimeStep, lambda Time, Velocity, Displacement: self.Acceleration(Displacement, Velocity)) 

def SimulateSpring(Mass, SpringConstant, DampingRatio, TimeStep):
    load = Load(Mass)
    CriticalDampingCoefficient = 2 * math.sqrt(SpringConstant * Mass)
    TimeConstant = 1 / math.sqrt(SpringConstant / Mass)
    Duration = 20 * TimeConstant
    DampingCoefficient = DampingRatio * CriticalDampingCoefficient
    springEuler = EulerSpring(0, SpringConstant, load, DampingCoefficient)
    springRungeKutta = RungeKuttaSpring(0, SpringConstant, load, DampingCoefficient)
    DataRungeKutta = []
    DataEuler = []
    Timestamps = []
    time = 0
    while time < Duration:
        springEuler.Update(TimeStep)
        springRungeKutta.Update(TimeStep)
        time += TimeStep
        Timestamps.append(time)    
        DataEuler.append(springEuler.Displacement.value)    
        DataRungeKutta.append(springRungeKutta.Displacement.value)    

    plt.plot(Timestamps, DataEuler)
    plt.plot(Timestamps, DataRungeKutta)
    plt.show()

def main(argv):
    Mass = float(argv[1])
    SpringConstant = float(argv[2])
    DampingCoefficient = float(argv[3])
    TimeStep = float(argv[4])
    SimulateSpring(Mass, SpringConstant, DampingCoefficient, TimeStep)

if __name__ == "__main__":
    main(sys.argv)

