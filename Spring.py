import math
import matplotlib.pyplot as plt
import sys

TimeStep = 0.001
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
    def __init__(self, InitialDisplacement, SpringConstant, Load):
        self.Displacement = PhysicalQuantity(InitialDisplacement, 0)
        self.SpringConstant = SpringConstant
        self.Load = Load 
        self.Force = self.Displacement.value * self.SpringConstant
        self.DampingCoefficient = 0 

    def Update(self, TimeStep):
        RK4SecondOrder(self.Displacement, TimeStep, lambda Time, Displacement, Velocity: self.Acceleration(Displacement, Velocity)) 

    def Acceleration(self, Displacement, Velocity):
        return (self.Load.Weight - self.DampingCoefficient * Velocity - self.SpringConstant * Displacement) / self.Load.Mass

def RK4SecondOrder(value, step, function, time = 0):
    derivativeCoefficients = []
    valueCoefficients = []

    derivativeCoefficients.append(step * function(time, value.value, value.derivative))
    valueCoefficients.append(step * value.derivative)

    for i in range(2):
        derivativeCoefficients.append(step * function(time + step / 2, value.derivative + derivativeCoefficients[-1] / 2, value.value + valueCoefficients[-1] / 2)) 
        valueCoefficients.append(step * (value.derivative + derivativeCoefficients[-1] / 2))

    derivativeCoefficients.append(step * function(time + step, value.derivative + derivativeCoefficients[-1], value.value + valueCoefficients[-1]))
    valueCoefficients.append(step * (value.derivative + derivativeCoefficients[-1]))

    value.derivative += (derivativeCoefficients[0] + 2 * derivativeCoefficients[1] + 2 * derivativeCoefficients[2] + derivativeCoefficients[3]) / 6
    value.value += (valueCoefficients[0] + 2 * valueCoefficients[1] + 2 * valueCoefficients[2] + valueCoefficients[3]) / 6
 
def SimulateSpring(Mass, SpringConstant, DampingCoefficient, Duration):
    load = Load(LoadKg)
    spring = Spring(0, SpringConstant, load)
    Data = []
    Timestamps = []
    time = 0
    while time < Duration:
        spring.Update(TimeStep)
        time += TimeStep
        Timestamps.append(time)    
        Data.append(spring.Displacement.value)    

    plt.plot(Timestamps, Data)
    plt.show()

def main(argv):
    Mass = float(argv[1])
    SpringConstant = float(argv[2])
    DampingCoefficient = float(argv[3])
    Duration = float(argv[4])
    SimulateSpring(Mass, SpringConstant, DampingCoefficient, Duration)

if __name__ == "__main__":
    main(sys.argv)

