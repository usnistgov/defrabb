import os
import sys  
from pybedtools import BedTool  
import pandas as pd  

def calculate_total_bp(bed):  
    """Calculate the total base pairs in a BED file."""  
    return sum(len(feature) for feature in bed)  
  
def main(dip_bed_path, summary_table_path, exclusion_paths):  
    # Read the dip.bed file  
    dip_bed = BedTool(dip_bed_path)  
    total_dip_bp = calculate_total_bp(dip_bed)  
  
    # Initialize a list to store summary info  
    summary_data = []  
  
    ## Adding initial dip size 
        # Add summary for remaining bases  
    summary_data.append({  
        'exclusion': 'dip.bed',  
        'excluded_bp': total_dip_bp,  
        'excluded_pct': 0,  
    }) 

    # Create an empty BedTool object for accumulating intersections  
    accumulated_intersections = None 

    # Use BedTool's multi_intersect to find intersections with all exclusions  
    if exclusion_paths:
        # Assuming you have a directory path where you want to save the files  
        output_dir = "intersections"  
        if not os.path.exists(output_dir):  
            os.makedirs(output_dir)  
            
        # Process each exclusion individually for their overlap with dip.bed  
        for exclusion_path in exclusion_paths:   
            # Calculate individual intersections  
            exclusion_bed = BedTool(exclusion_path)  
            intersection = dip_bed.intersect(exclusion_bed)

            intersection_bp = calculate_total_bp(intersection)  
            intersection_pct = (intersection_bp / total_dip_bp) * 100

            # Extract the exclusion filename for summary  
            exclusion_filename = exclusion_path.split('/')[-1]  

            # Save the intersection to a file before accumulating  
            intersection_filename = os.path.join(output_dir, exclusion_filename)  
            intersection.saveas(intersection_filename)  
  
            summary_data.append({  
                'exclusion': exclusion_filename,  
                'excluded_bp': intersection_bp,  
                'excluded_pct': intersection_pct,  
            })  

            # Update accumulated_intersections with the saved file  
            if accumulated_intersections is None:  
                accumulated_intersections = BedTool(intersection_filename)  
            else:  
                accumulated_intersections = accumulated_intersections.cat(intersection_filename, postmerge=True)

        # Merge accumulated intersections to consolidate overlaps
        total_excluded_bp = calculate_total_bp(accumulated_intersections)
        total_excluded_pct = (total_excluded_bp/ total_dip_bp) * 100
        remaining_bp = total_dip_bp - total_excluded_bp  
        remaining_pct = (remaining_bp / total_dip_bp) * 100
    else:  
        # If there were no exclusions  
        remaining_bp = total_dip_bp  
        remaining_pct = 100
  
    # Add summary for remaining bases  
    summary_data.append({  
        'exclusion': 'benchmark_regions',  
        'excluded_bp': remaining_bp,  
        'excluded_pct': remaining_pct,  
    })  
  
    # Add summary for cumulative intersections  
    summary_data.append({  
        'exclusion': 'total_excluded',  
        'excluded_bp': total_excluded_bp,  
        'excluded_pct': total_excluded_pct,  
    })  
  
    # Convert summary data to DataFrame and save  
    summary_df = pd.DataFrame(summary_data)  
    summary_df.to_csv(summary_table_path, index=False)  
  
if __name__ == "__main__":
    print(sys.argv[1])
    dip_bed_path = sys.argv[1]  
    summary_table_path = sys.argv[2]  
    exclusion_paths = sys.argv[3:]  
    main(dip_bed_path, summary_table_path, exclusion_paths)  

