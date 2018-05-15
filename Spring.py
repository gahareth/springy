import math
import matplotlib.pyplot as plt
import sys
from Core import AccelerationDueToGravity
from Core import RK4SecondOrder
from Core import Euler 
from Core import PhysicalQuantity
from Load import Load

class ComplexNumber():
    def __init__(self, real, imaginary):
        self.real = real
        self.imaginary = imaginary

class Spring:
    def __init__(self, InitialDisplacement, SpringConstant, Load, DampingCoefficient):
        self.EquilibriumDisplacement = Load.Weight / SpringConstant  
        self.InitialDisplacement = InitialDisplacement
        self.InitialVelocity = InitialDisplacement - self.EquilibriumDisplacement
        self.Displacement = PhysicalQuantity(InitialDisplacement, 0)
        self.SpringConstant = SpringConstant
        self.Load = Load 
        self.DampingCoefficient = DampingCoefficient
        self.CalculatePoles()

    def AnalyticalSolution(self, time):
        return self.EquilibriumDisplacement - (self.EquilibriumDisplacement * math.exp(self.Poles[0].real * time) * math.cos(self.Poles[0].imaginary * time) + (self.InitialVelocity / self.Poles[0].real) * math.exp(self.Poles[1].real * time) * math.sin(self.Poles[1].imaginary * time))
        
    def CalculatePoles(self):
        discriminant = self.DampingCoefficient ** 2 - 4 * self.Load.Mass * self.SpringConstant
        if discriminant <= 0:
            poleRealPart = -self.DampingCoefficient / (2 * self.Load.Mass)
            poleImaginaryPart = math.sqrt(-discriminant) / (2 * self.Load.Mass)
            self.Poles = [ComplexNumber(poleRealPart, poleImaginaryPart), ComplexNumber(poleRealPart, -poleImaginaryPart)]
        else:
            pole1 = (-self.DampingCoefficient + math.sqrt(discriminant)) / (2 * self.Load.Mass)
            pole2 = (-self.DampingCoefficient - math.sqrt(discriminant)) / (2 * self.Load.Mass)
            self.Poles = [ComplexNumber(pole1, 0), ComplexNumber(pole2, 0)]
        
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
    DataAnalytical = []
    DataRungeKutta = []
    DataEuler = []
    Timestamps = []
    time = 0
    while time < Duration:
        springEuler.Update(TimeStep)
        springRungeKutta.Update(TimeStep)
        time += TimeStep
        Timestamps.append(time)    
        DataAnalytical.append(springEuler.AnalyticalSolution(time))
        DataEuler.append(springEuler.Displacement.value)    
        DataRungeKutta.append(springRungeKutta.Displacement.value)    

    plt.plot(Timestamps, DataAnalytical, label='Actual')
    plt.plot(Timestamps, DataEuler, label='Euler')
    plt.plot(Timestamps, DataRungeKutta, label='Runge-Kutta')
    plt.xlabel('time (s)')
    plt.ylabel('displacement (m)')
    plt.legend()
    plt.show()

def main(argv):
    Mass = float(argv[1])
    SpringConstant = float(argv[2])
    DampingCoefficient = float(argv[3])
    TimeStep = float(argv[4])
    SimulateSpring(Mass, SpringConstant, DampingCoefficient, TimeStep)

if __name__ == "__main__":
    main(sys.argv)

