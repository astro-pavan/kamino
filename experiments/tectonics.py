import argparse
import itertools
import os
from kamino.planet2 import Planet
import kamino.planet2 as p2
from kamino.constants import M_EARTH, R_EARTH

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Sweep parameters for planet tectonics.")
    parser.add_argument('--reverse-weathering', action='store_true', 
                        help="Enable reverse weathering.")
    parser.add_argument('--crust-carbonate-fraction', type=float, default=0.0, 
                        help="Set the carbonate crust fraction (default: 0.0).")
    args = parser.parse_args()

    # Override the output path from planet2
    output_path = '/home/pavan/PhD/kamino_experiments'
    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path

    # Create the output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    # Standard parameters used for the Planet initialization
    BACKGROUND_PRESSURE = 1e5   # Pa (1 bar)
    OCEAN_DEPTH = 3000          # m

    # Parameter sweeps
    instellation = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
    outgassing = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    crust_production_rate = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]

    # Iterate through all combinations of the parameters
    for s, o, c in itertools.product(instellation, outgassing, crust_production_rate):
        
        # Construct a unique name for each simulation
        run_name = f'planet_s_{s}_out_{o}_crust_{c}'
        if args.reverse_weathering:
            run_name += '_rw'
        if args.crust_carbonate_fraction > 0.0:
            run_name += f'_carb_{args.crust_carbonate_fraction}'
        
        print(f"--- Starting: {run_name} ---")
        print(f"Instellation: {s}, Outgassing: {o}, Crust Production: {c}")
        print(f"Reverse Weathering: {args.reverse_weathering}, Crust Carbonate Content: {args.crust_carbonate_fraction}")
        
        try:
            # Initialize the planet with the runtime arguments
            p = Planet(
                mass=M_EARTH,
                radius=R_EARTH,
                background_pressure=BACKGROUND_PRESSURE,
                instellation=s,
                crust_production_rate=c,
                outgassing=o,
                ocean_depth=OCEAN_DEPTH,
                crust_carbonate_content=args.crust_carbonate_fraction,
                reverse_weathering=args.reverse_weathering,
                name=run_name
            )
            
            # Run the simulation
            p.time_evolve()
            
        except Exception as e:
            print(f"Simulation failed for {run_name} with error: {e}")
            
        print(f"--- Finished: {run_name} ---\n")


if __name__ == "__main__":
    main()