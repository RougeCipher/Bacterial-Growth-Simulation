import numpy as np
from .advanced_diffusion import adi_diffusion_2d_variableD

class SpatialEnvironment:
    """
    2D environment for ADI-based approach:
      - nutrient_grid, waste_grid
      - boundary conditions
      - optional variable diffusion
      - time-dependent environment states (T, pH, etc.)
    """
    def __init__(self, env_config):
        self.nx = env_config.grid_size_x
        self.ny = env_config.grid_size_y
        self.dx = env_config.dx
        self.dy = env_config.dy
        self.bc_type = env_config.boundary_condition_type
        self.bc_value= env_config.boundary_value

        self.temperature = env_config.initial_temperature
        self.pH = env_config.initial_pH
        self.oxygen_level= env_config.oxygen_level
        self.moisture_level= env_config.moisture_level
        self.nutrients_dict= env_config.nutrients
        self.max_waste= env_config.max_waste

        if env_config.variable_diffusion and env_config.diffusion_map is not None:
            self.D = np.array(env_config.diffusion_map, dtype=np.float64)
            if self.D.shape!=(self.ny, self.nx):
                raise ValueError("diffusion_map shape mismatch with grid_size_x,grid_size_y.")
        else:
            self.D = float(env_config.diffusion_coefficient)

        self.nutrient_grid= np.zeros((self.ny, self.nx), dtype=np.float64)
        key, val = next(iter(self.nutrients_dict.items()))
        self.nutrient_grid[:,:] = val
        self.waste_grid= np.zeros((self.ny, self.nx), dtype=np.float64)

    def update_env_state(self, t, dt):
        """
        For demonstration, sinusoidal variations in T and pH.
        """
        self.temperature += 0.5*np.sin(2*np.pi*t/24)
        self.pH += 0.05*np.sin(2*np.pi*t/24)

    def diffuse_nutrients(self, dt, consumption):
        """
        ADI approach to diffuse the nutrient field with a reaction = -consumption
        """
        self.nutrient_grid = adi_diffusion_2d_variableD(
            conc=self.nutrient_grid,
            D=self.D,
            dt=dt,
            dx=self.dx,
            dy=self.dy,
            boundary_type=self.bc_type,
            boundary_value=self.bc_value,
            reaction=-consumption
        )

    def accumulate_waste(self, dt, production):
        """
        Add waste. Then clamp at max_waste.
        """
        self.waste_grid += production*dt
        np.clip(self.waste_grid, 0, self.max_waste, out=self.waste_grid)







