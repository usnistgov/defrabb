import argparse 
import os
from typing import List
from pybedtools import BedTool  
import pandas as pd
import logging 

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')  

import argparse  
import os  
from typing import List  
from pybedtools import BedTool  
import pandas as pd  
import logging  
  
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')  
  
  
def calculate_total_bp(bed: BedTool) -> int:  
    """Calculate the total base pairs in a BED file."""  
    return sum(len(feature) for feature in bed)  
  
  
def process_exclusion(exclusion_path: str, dip_bed: BedTool, output_dir: str) -> (str, int, float):  
    """Process individual exclusion file to find and save intersection."""  
    exclusion_bed = BedTool(exclusion_path)  
    intersection = dip_bed.intersect(exclusion_bed)  
    intersection_bp = calculate_total_bp(intersection)  
    exclusion_filename = os.path.basename(exclusion_path)  
    intersection_filename = os.path.join(output_dir, exclusion_filename)  
    intersection.saveas(intersection_filename)  
    logging.info(f"Processed exclusion: {exclusion_filename}")  
    return intersection_filename, intersection_bp  
  
  
def accumulate_intersections(exclusion_filenames: List[str]) -> BedTool:  
    """Accumulate intersections from individual exclusion files."""  
    accumulated_intersections = BedTool(exclusion_filenames[0])  
    for filename in exclusion_filenames[1:]:  
        accumulated_intersections = accumulated_intersections.cat(filename, postmerge=True)  
    return accumulated_intersections  
  
  
def generate_summary_data(dip_bed_path: str, exclusion_paths: List[str], intersection_dir: str) -> List[dict]:  
    """Generate summary data for exclusions."""  
    summary_data = []  
    dip_bed = BedTool(dip_bed_path)  
    total_dip_bp = calculate_total_bp(dip_bed)  
    summary_data.append({'genomic_region': 'initial', 'bp': total_dip_bp, 'pct_of_initial': 0})  
  
    if exclusion_paths:   
        processed_exclusions = [process_exclusion(path, dip_bed, intersection_dir) for path in exclusion_paths]  
        exclusion_filenames, intersections_bp = zip(*processed_exclusions)  
        intersection_pct = [(bp / total_dip_bp) * 100 for bp in intersections_bp]  
  
        for path, bp, pct in zip(exclusion_paths, intersections_bp, intersection_pct):  
            summary_data.append({  
                'genomic_region': os.path.basename(path),  
                'bp': bp,  
                'pct_of_initial': pct,  
            })  
  
        accumulated_intersections = accumulate_intersections(exclusion_filenames)  
        total_excluded_bp = calculate_total_bp(accumulated_intersections)  
        total_excluded_pct = (total_excluded_bp / total_dip_bp) * 100  
        remaining_bp = total_dip_bp - total_excluded_bp  
        remaining_pct = (remaining_bp / total_dip_bp) * 100  
  
        summary_data.append({'genomic_region': 'benchmark_regions', 'bp': remaining_bp, 'pct_of_initial': remaining_pct})  
        summary_data.append({'genomic_region': 'total_excluded', 'bp': total_excluded_bp, 'pct_of_initial': total_excluded_pct})  
    else:  
        logging.info("No exclusions provided.")  
        remaining_bp = total_dip_bp  
        remaining_pct = 100  
        summary_data.append({'genomic_region': 'benchmark_regions', 'bp': remaining_bp, 'pct_of_initial': remaining_pct})  
  
    return summary_data  
  
  
def save_summary_table(summary_data: List[dict], summary_table_path: str):  
    """Save summary data to a CSV file."""  
    summary_df = pd.DataFrame(summary_data)  
    summary_df.to_csv(summary_table_path, index=False)  
    logging.info(f"Summary table saved to {summary_table_path}")  
  
  
def parse_arguments() -> argparse.Namespace:  
    parser = argparse.ArgumentParser(description='Calculate the excluded base pairs from a dip.bed file using exclusion bed files.')  
    parser.add_argument('dip_bed_path', type=str, help='Path to the dip.bed file.')  
    parser.add_argument('summary_table_path', type=str, help='Path where the summary table will be saved.') 
    parser.add_argument('intersection_dir', type=str, help='Directory to write dip.bed and exclusion intersection bed files.')   
    parser.add_argument('exclusion_paths', type=str, nargs='+', help='Paths to the exclusion bed files.')  
    return parser.parse_args()  
  
  
def main(dip_bed_path: str, summary_table_path: str, exclusion_paths: List[str],intersection_dir: str):  
    logging.info("Starting script.")  
    os.makedirs(intersection_dir, exist_ok=True) 
    summary_data = generate_summary_data(dip_bed_path, exclusion_paths, intersection_dir)  
    save_summary_table(summary_data, summary_table_path)  
    logging.info("Script finished.")  
  
  
if __name__ == "__main__":  
    args = parse_arguments()  
    main(args.dip_bed_path, args.summary_table_path, args.exclusion_paths, args.intersection_dir)  
  


