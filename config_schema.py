from pydantic import BaseModel, Field, validator
from typing import Dict, Optional

class AntibioticConfig(BaseModel):
    name: str = Field("Ciprofloxacin", description="Antibiotic name.")
    enabled: bool = Field(True, description="If True, antibiotic is used.")
    initial_efficiency: float = Field(0.95, ge=0, le=1, description="Initial antibiotic efficiency [0..1].")
    decay_rate: float = Field(0.02, ge=0, description="Exponential decay rate of antibiotic efficiency.")

class BacteriumConfig(BaseModel):
    species_name: str = Field("Escherichia coli", description="Name of the bacterium.")
    base_growth_rate: float = Field(0.25, ge=0, description="Base growth rate (hr^-1) under optimal conditions.")
    optimal_temp: float = Field(37.0, description="Optimal temperature for the bacterium (°C).")
    optimal_pH: float = Field(7.0, description="Optimal pH for the bacterium.")

class EnvironmentConfig(BaseModel):
    grid_size_x: int = Field(50, description="Number of cells in X dimension.")
    grid_size_y: int = Field(50, description="Number of cells in Y dimension.")
    dx: float = Field(1.0, description="Spatial step in X dimension (e.g., 1 mm).")
    dy: float = Field(1.0, description="Spatial step in Y dimension (e.g., 1 mm).")
    diffusion_coefficient: float = Field(0.05, ge=0, description="Diffusion coefficient for nutrients (in whatever units).")
    initial_temperature: float = 37.0
    initial_pH: float = 7.0
    oxygen_level: float = 1.0
    moisture_level: float = 1.0
    nutrients: Dict[str, float] = Field(default_factory=lambda: {"glucose": 1.0})
    max_waste: float = 10.0

class SimulationConfig(BaseModel):
    time_steps: int = Field(24, ge=1, description="Number of time steps to simulate.")
    dt: float = Field(1.0, gt=0, description="Time step in hours.")
    initial_population: float = Field(1000, ge=0, description="Initial total population (non-resistant + resistant).")
    resistant_population: float = Field(10, ge=0, description="Initial resistant subpopulation.")
    carrying_capacity: float = Field(1e9, ge=1, description="Logistic carrying capacity.")

class NCBIConfig(BaseModel):
    api_key: str = Field("", description="NCBI API key.")
    email: str = Field("", description="User email associated with NCBI API.")

class RootConfig(BaseModel):
    antibiotic: AntibioticConfig
    bacterium: BacteriumConfig
    environment: EnvironmentConfig
    simulation: SimulationConfig
    ncbi: NCBIConfig

    @validator("simulation")
    def check_population(cls, v):
        if v.resistant_population > v.initial_population:
            raise ValueError("Resistant population cannot exceed total initial population.")
        return v
