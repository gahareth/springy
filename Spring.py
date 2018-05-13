import math
import matplotlib.pyplot as plt
import sys

LoadKg = 10

AccelerationDueToGravity = 9.8

class Load:
    def __init__(self, Mass):
        self.Mass = Mass
        self.Weight = Mass * AccelerationDueToGravity

class PhysicalQuantity:
    def __init__(self, InitialValue, InitialDerivative):
        self.value = InitialValue
        self.derivative = InitialDerivative

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

def Euler(value, step, function, time = 0):
    initialValue = value.value
    value.value += step * value.derivative
    value.derivative += step * function(time, value.derivative, initialValue)

def RK4SecondOrder(value, step, function, time = 0):
    derivativeCoefficients = []
    valueCoefficients = []

    derivativeCoefficients.append(step * function(time, value.derivative, value.value))
    valueCoefficients.append(step * value.derivative)

    for i in range(2):
        derivativeCoefficients.append(step * function(time + step / 2, value.derivative + derivativeCoefficients[-1] / 2, value.value + valueCoefficients[-1] / 2)) 
        valueCoefficients.append(step * (value.derivative + derivativeCoefficients[-1] / 2))

    derivativeCoefficients.append(step * function(time + step, value.derivative + derivativeCoefficients[-1], value.value + valueCoefficients[-1]))
    valueCoefficients.append(step * (value.derivative + derivativeCoefficients[-1]))

    value.derivative += (derivativeCoefficients[0] + 2 * derivativeCoefficients[1] + 2 * derivativeCoefficients[2] + derivativeCoefficients[3]) / 6
    value.value += (valueCoefficients[0] + 2 * valueCoefficients[1] + 2 * valueCoefficients[2] + valueCoefficients[3]) / 6
 
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

