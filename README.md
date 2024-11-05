# Bacterial Growth Simulation

This project is a Python simulation of bacterial growth under various environmental and antibiotic conditions. It models the growth of bacterial populations over time, taking into account factors like temperature, pH, nutrient availability, and the impact of antibiotics on resistant and non-resistant populations.

## Features

- Models bacterial population dynamics with antibiotic resistance.
- Includes environmental factors like temperature, pH, oxygen, and moisture levels.
- Simulates the depletion of nutrients and accumulation of waste.
- Provides graphical visualization of:
  - Total, non-resistant, and resistant bacterial populations over time.
  - Antibiotic efficiency decay over time.
  - Environmental conditions (temperature and pH) over time.
  - Waste accumulation and nutrient depletion over time.

## Requirements

- Python 3.7 or higher
- Libraries:
  - `matplotlib`
  - `numpy`
  - `tkinter` (for user inputs, included with Python)

To install the necessary packages, run:
```bash
pip install matplotlib numpy

## Usage

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/Bacterial-Growth-Simulation.git
   cd Bacterial-Growth-Simulation
   ```

2. **Run the simulation**:
   Execute the main script using Python:
   ```bash
   python main.py
   ```

3. **Enter simulation parameters**:
   - The program will prompt you to enter parameters like bacterial species, antibiotic type, nutrient environment, temperature, pH, oxygen level, and more.
   - Example inputs:
     - Bacterial species: `Escherichia coli`
     - Antibiotic: `Ciprofloxacin`
     - Nutrient environment: `Rich`
     - Initial population size: `1000`
     - Timeframe for simulation: `24` hours

4. **View Results**:
   - The simulation outputs a series of plots showing bacterial growth, antibiotic decay, population diversity, environmental conditions, waste accumulation, and nutrient depletion.
   - A summary of key metrics, such as phase durations and peak population, is also displayed.

## Simulation Summary Example

An example of the output summary:

```
---------------------------------------
Simulation Summary
---------------------------------------
Lag Phase Duration: 7.0 hours
Peak Population: 2.9 x 10^3 CFU
Phase Durations:
  - Lag: 7.0 hours
  - Exponential: 10.0 hours
  - Stationary: 0.0 hours
  - Death: 7.0 hours
Environmental Summary:
  - Temperature variation: ±5°C around 37°C
  - pH variation: ±0.3 around pH 7.0
  - Nutrient Depletion Rate: 0.00001
  - Waste Accumulation: Growth-inhibiting after 24.0 hours
```

## Graphical Output

The simulation generates the following plots:
- **Bacterial Growth Over Time**: Total, non-resistant, and resistant populations.
- **Antibiotic Efficiency Decay**: Displays the decay of antibiotic effectiveness.
- **Population Diversity**: Shows non-resistant vs. resistant population changes.
- **Environmental Conditions**: Temperature and pH changes over time.
- **Waste Accumulation**: Tracks waste levels as the simulation progresses.
- **Nutrient Depletion**: Shows how nutrient levels decrease over time.

## License

This project is licensed under the MIT License.

## Acknowledgements

This project was inspired by microbial population dynamics and the study of antibiotic resistance in controlled environments.
```

