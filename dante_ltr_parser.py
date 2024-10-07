import argparse
import subprocess
import os
import traceback
from python_scripts.generate_csv import TEAnalyzer
from python_scripts.generate_report import generate_report

# args: genome fasta, dante_LTR GFF3, alternativelly s. substitution rate
def args_from_parser():
    parser = argparse.ArgumentParser(
        description='''Script performs analysis of individual TEs in DANTE_LTR GFF3 outputs''')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        "-fasta", "--genome_fasta", type=str, required=True,
        help='Genome FASTA file'
        )
    requiredNamed.add_argument(
        "-gff", "--dante_ltr_gff", type=str, required=True,
        help='GFF3 output file from DANTE_LTR'
        )
    requiredNamed.add_argument(
        "-spec", "--species", type=str, required=True,
        help='Latin name of studied species \'Arabidopsis thaliana\''
        )
    parser.add_argument(
        "-ssr", "--synonymous_substitution_rate", type=float, default=1.5e-8,
        help='If not specified SSR for Arabidopsis: 1.5e-8 will be used'
        )
    parser.add_argument(
        "-t", "--threads", type=int, default=os.cpu_count() / 2,
        help=f'Number of CPUs for LTR sequences analysis default is {int(os.cpu_count() / 2)} CPU for your PC (max of CPU count / 2)'
        )
    parser.add_argument(
        "-out_path", "--output_path", type=str, default='/home/',
        help='Path where the output should be stored'
        )
    return parser.parse_args()

import subprocess

def main():
    args = args_from_parser()
    
    try:
        # Process GFF3 and generate CSV
        analyzer = TEAnalyzer(args.dante_ltr_gff, args.genome_fasta, args.synonymous_substitution_rate, args.output_path, args.threads)
        out_csv_path = analyzer.csv_generator()
        # out_csv_path = "/home/pavel/Documents/Work/42_LH_Humulus_satellites/Akagi_Hj_Hl_haploidGenomeAssemeblies_Sep_2024/10-12_haploid/DANTE_LTR_summary/Hop_10_12hap_LTR_retrotransposons_annotation_chrOnly_TE_characteristics.csv"
        
        # Generate plots using R script
        cmd = f"Rscript r_scripts/Rplots_from_csv.R {out_csv_path}"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"RSCRIPT process failed with return code {process.returncode}")
            print(f"Error: {stderr.decode('utf-8')}")
        else:
            print(f"RSCRIPT process completed successfully. Output: {stdout.decode('utf-8')}")
        
        # Generate report
        generate_report(args.output_path, args.species)
    
    except Exception as e:
        print(f"An error occurred: {e}")
    # Print detailed traceback information
        traceback.print_exc()

if __name__ == "__main__":
    main()