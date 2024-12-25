import numpy as np

class Antibiotic:
    """
    Basic antibiotic model with exponential decay over time and 
    environment-based efficiency penalties.
    """

    def __init__(self, name, initial_efficiency, decay_rate):
        self.name = name
        self.initial_efficiency = initial_efficiency
        self.decay_rate = decay_rate

    def efficiency(self, t, temperature, pH):
        """
        Returns the current antibiotic efficiency (0..1) at time t, 
        factoring in exponential decay and environment-based penalties.
        """
        eff = self.initial_efficiency * np.exp(-self.decay_rate * t)

        # Temperature penalty if too cold/hot
        if temperature < 15 or temperature > 45:
            eff *= 0.7

        # pH penalty if too acidic/basic
        if pH < 5 or pH > 9:
            eff *= 0.8

        return float(np.clip(eff, 0, 1))
