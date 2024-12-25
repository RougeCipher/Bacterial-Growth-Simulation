import numpy as np

class Antibiotic:
    def __init__(self, name, initial_efficiency, decay_rate):
        self.name = name
        self.initial_efficiency = initial_efficiency
        self.decay_rate = decay_rate

    def efficiency(self, t, temperature, pH):
        eff = self.initial_efficiency * np.exp(-self.decay_rate * t)
        if temperature < 15 or temperature > 45:
            eff *= 0.7
        if pH < 5 or pH > 9:
            eff *= 0.8
        return float(np.clip(eff, 0, 1))
