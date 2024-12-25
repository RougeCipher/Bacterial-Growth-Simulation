import numpy as np

class BacteriaGrowth2D:
    """
    Bacterial population model with logistic growth + antibiotic.
    Each cell in 2D has non-resistant + resistant populations.
    """

    def __init__(self, config, env_config):
        self.species_name = config.species_name
        self.base_growth_rate = config.base_growth_rate
        self.optimal_temp = config.optimal_temp
        self.optimal_pH = config.optimal_pH

        # For 2D: get environment dimension
        self.nx = env_config.grid_size_x
        self.ny = env_config.grid_size_y
        self.carrying_capacity = 1e9

        total_pop = config.simulation.initial_population
        resistant_pop = config.simulation.resistant_population
        fraction_resistant = resistant_pop / total_pop if total_pop > 0 else 0

        uniform_density = total_pop / (self.nx * self.ny)
        uniform_res_density = uniform_density * fraction_resistant
        uniform_nonres_density = uniform_density * (1 - fraction_resistant)

        self.non_resistant_grid = np.full((self.ny, self.nx), uniform_nonres_density, dtype=np.float64)
        self.resistant_grid = np.full((self.ny, self.nx), uniform_res_density, dtype=np.float64)

    def apply_environment_factors_cell(self, env, y, x):
        """
        Combine temperature, pH, oxygen, moisture, local nutrient, and local waste 
        to adjust the base growth rate in a single cell (y,x).
        """
        gr = self.base_growth_rate

        # Temperature penalty
        if abs(env.temperature - self.optimal_temp) > 5:
            gr *= 0.5
        elif abs(env.temperature - self.optimal_temp) > 2:
            gr *= 0.8

        # pH penalty
        if abs(env.pH - self.optimal_pH) > 1:
            gr *= 0.7

        # Oxygen
        gr *= np.clip(env.oxygen_level, 0, 1)

        # Moisture
        gr *= np.clip(env.moisture_level, 0, 1)

        # Nutrients
        nutr = env.nutrient_grid[y, x]
        gr *= nutr

        # Waste penalty
        waste_cell = env.waste_grid[y, x]
        factor = 1 - (waste_cell / env.max_waste)
        gr *= max(factor, 0.0)

        return gr

    def update_step(self, dt, env, antibiotic_eff):
        """
        Performs one time-step update using a forward Euler approach in each cell.
        Returns a 2D array 'consumption_grid' indicating nutrient consumption 
        (based on the amount of bacterial growth).
        """
        consumption_grid = np.zeros_like(self.non_resistant_grid)

        for y in range(self.ny):
            for x in range(self.nx):
                nr = self.non_resistant_grid[y, x]
                r = self.resistant_grid[y, x]
                total = nr + r

                # If nearly zero, skip
                if total < 1e-12:
                    continue

                env_gr = self.apply_environment_factors_cell(env, y, x)

                # Logistic factor
                logistic_term = (1 - total / self.carrying_capacity)

                # Non-resistant growth (affected by antibiotic)
                nr_growth_rate = env_gr * (1 - antibiotic_eff) * logistic_term

                # Resistant growth (less/no impact from antibiotic)
                r_growth_rate = env_gr * logistic_term

                d_nr = dt * nr * nr_growth_rate
                d_r = dt * r * r_growth_rate

                # Mutation from non-resistant to resistant
                mutation_rate = 0.0001 + 0.001 * antibiotic_eff
                mutated = nr * mutation_rate * dt

                nr_new = nr + d_nr - mutated
                r_new = r + d_r + mutated

                # No negatives
                if nr_new < 0:
                    nr_new = 0
                if r_new < 0:
                    r_new = 0

                self.non_resistant_grid[y, x] = nr_new
                self.resistant_grid[y, x] = r_new

                # Nutrient consumption is proportionate to bacterial growth
                consumption_grid[y, x] = d_nr + d_r

        return consumption_grid
