from pydantic import BaseModel, Field, validator
from typing import Dict, Optional, List

class AntibioticConfig(BaseModel):
    name: str = "Ciprofloxacin"
    enabled: bool = True
    initial_efficiency: float = Field(0.95, ge=0, le=1)
    decay_rate: float = Field(0.02, ge=0)

class BacteriumConfig(BaseModel):
    species_name: str = "Escherichia coli"
    base_growth_rate: float = Field(0.25, ge=0)
    optimal_temp: float = 37.0
    optimal_pH: float = 7.0

class EnvironmentConfig(BaseModel):
    grid_size_x: int = 50
    grid_size_y: int = 50
    dx: float = 1.0
    dy: float = 1.0
    diffusion_coefficient: float = 0.05
    variable_diffusion: bool = False
    diffusion_map: Optional[List[List[float]]] = None
    boundary_condition_type: str = "dirichlet"
    boundary_value: float = 0.0
    initial_temperature: float = 37.0
    initial_pH: float = 7.0
    oxygen_level: float = 1.0
    moisture_level: float = 1.0
    nutrients: Dict[str, float] = {"glucose": 1.0}
    max_waste: float = 10.0

class SimulationConfig(BaseModel):
    time_steps: int = 24
    dt: float = 1.0
    initial_population: float = 1000
    resistant_population: float = 10
    carrying_capacity: float = 1e9

class NCBIConfig(BaseModel):
    api_key: str = ""
    email: str = ""

class RootConfig(BaseModel):
    antibiotic: AntibioticConfig
    bacterium: BacteriumConfig
    environment: EnvironmentConfig
    simulation: SimulationConfig
    ncbi: NCBIConfig

    @validator("simulation")
    def check_populations(cls, v):
        if v.resistant_population > v.initial_population:
            raise ValueError("resistant_population cannot exceed initial_population")
        return v
