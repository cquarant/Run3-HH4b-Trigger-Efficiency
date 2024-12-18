import argparse
from pathlib import Path
import uproot

years = ["2022", "2023"]  # TODO: add more years if needed

def calculate_kfactor(year):
    """Calculate the k-factor for a given year using only FatJet1_pt."""
    path_data = Path(f"hists/Histograms_{year}_data.root")
    MC_regions = ["QCD", "TTbar", "VJ", "VV"]
    paths_MC = {
        region: Path(f"hists/Histograms_{year}_MC_{region}.root")
        for region in MC_regions
    }

    # First, load data histogram
    try:
        data_file = uproot.open(path_data)
    except FileNotFoundError:
        raise FileNotFoundError(f"Data file for year {year} not found")

    data_hist = data_file["FatJet1_pt;1"]
    data_integral = data_hist.values().sum()

    # Combine MC histograms
    combined_mc = None
    for region in MC_regions:
        mc_file = uproot.open(paths_MC[region])
        hist = mc_file["FatJet1_pt;1"]
        
        if combined_mc is None:
            combined_mc = hist.values()
        else:
            combined_mc += hist.values()

    # Calculate data/MC ratio
    mc_integral = combined_mc.sum()
    return data_integral / mc_integral

def update_kfactor_file(year, kfactor):
    """Update or create kfact.h file."""
    script_dir = Path(__file__).parent
    kfactor_file = script_dir / "kfact.h"
    
    content = {}
    
    # Initialize with zeros
    for yr in years:
        content[yr] = 0.0
    
    # Read existing content if file exists
    if kfactor_file.exists():
        existing_content = kfactor_file.read_text()
        for yr in years:
            for line in existing_content.split('\n'):
                if f"#define KFACT_{yr}" in line:
                    try:
                        content[yr] = float(line.split()[-1])
                    except (ValueError, IndexError):
                        content[yr] = 0.0
    
    # Update the current year's k-factor
    content[year] = kfactor
    
    # Write all years to file with TODO comment
    header = "// TODO: Check if all k-factors are properly calculated\n"
    header += "// K-factors are calculated using FatJet1_pt distribution\n"
    header += "// A value of 0.000 means that year's k-factor hasn't been calculated yet\n\n"
    
    lines = [f"#define KFACT_{yr} {content[yr]:.3f}" for yr in years]
    kfactor_file.write_text(header + '\n'.join(lines) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Calculate k-factor for a specific year')
    parser.add_argument('--year', type=str, required=True, help='Year to process (e.g., "2022" or "2023")')
    args = parser.parse_args()
    
    try:
        kfactor = calculate_kfactor(args.year)
        update_kfactor_file(args.year, kfactor)
        print(f"Successfully calculated and updated k-factor for {args.year}: {kfactor:.3f}")
    except Exception as e:
        print(f"Error processing year {args.year}: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()