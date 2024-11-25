import json
import argparse
from pathlib import Path

def convert_json_to_txt(input_data):
    output_lines = []
    
    # Get sorted run numbers
    sorted_runs = sorted(input_data.keys(), key=int)
    
    # Process each run number
    for i, run_number in enumerate(sorted_runs):
        # Convert the ranges to the required format
        ranges = input_data[run_number]
        formatted_ranges = ",".join([f"{{{start}, {end}}}" for start, end in ranges])
        
        # Create the line in the required format
        # Add comma only if it's not the last line
        if i < len(sorted_runs) - 1:
            line = f"{{{run_number}, {{{formatted_ranges}}}}},"
        else:
            line = f"{{{run_number}, {{{formatted_ranges}}}}}"
        
        output_lines.append(line)
    
    # Join all lines with newlines and add an extra newline at the end
    return "\n".join(output_lines) + "\n"

def main():
    parser = argparse.ArgumentParser(description='Convert JSON input to specified text format')
    parser.add_argument('--input', required=True, type=Path, help='Input JSON file path')
    parser.add_argument('--output', required=True, type=Path, help='Output text file path')
    
    args = parser.parse_args()
    
    try:
        # Read input JSON file
        input_data = json.loads(args.input.read_text())
        
        # Convert to desired format
        output_text = convert_json_to_txt(input_data)
        
        # Ensure output directory exists and write output file
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output_text)
            
    except json.JSONDecodeError:
        print("Error: Invalid JSON file")
        exit(1)
    except IOError as e:
        print(f"Error: File operation failed - {str(e)}")
        exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()