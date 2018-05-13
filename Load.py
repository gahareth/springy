import math
import matplotlib.pyplot as plt
import sys
from Core import AccelerationDueToGravity
from Core import RK4SecondOrder
from Core import Euler 
from Core import PhysicalQuantity

class Load:
    def __init__(self, Mass):
        self.Mass = Mass
        self.Weight = Mass * AccelerationDueToGravity

    def Acceleration(self):
        return AccelerationDueToGravity

class RungeKuttaMass(Load):
    def __init__(self, Mass):
        Load.__init__(self, Mass)
        self.Displacement = PhysicalQuantity(0, 0)

    def Update(self, TimeStep):
        RK4SecondOrder(self.Displacement, TimeStep, lambda Time, Velocity, Displacement: self.Acceleration()) 

class EulerMass(Load):
    def __init__(self, Mass):
        Load.__init__(self, Mass)
        self.Displacement = PhysicalQuantity(0, 0)

    def Update(self, TimeStep):
        Euler(self.Displacement, TimeStep, lambda Time, Velocity, Displacement: self.Acceleration()) 

def SimulateMass(Mass, TimeStep, Duration):
    load = Load(Mass)
    DataRungeKutta = []
    DataEuler = []
    Timestamps = []
    time = 0
    massRungeKutta = RungeKuttaMass(Mass)
    massEuler = EulerMass(Mass)
    while time < Duration:
        massEuler.Update(TimeStep)
        massRungeKutta.Update(TimeStep)
        time += TimeStep
        Timestamps.append(time)    
        DataEuler.append(massEuler.Displacement.value)    
        DataRungeKutta.append(massRungeKutta.Displacement.value)    

    plt.plot(Timestamps, DataEuler)
    plt.plot(Timestamps, DataRungeKutta)
    plt.show()

def main(argv):
    Mass = float(argv[1])
    TimeStep = float(argv[2])
    Duration = float(argv[3])
    SimulateMass(Mass, TimeStep)

if __name__ == "__main__":
    main(sys.argv)
