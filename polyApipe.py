#!/usr/bin/env python3

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Count reads in polyA peaks. "+
     "Given a set of annotated bams files (cell barcode, UMI code and gene) "+
     "(1) Find where peaks of polyA reads occur. "+
     "(2) Count how many reads per sample fall in each polyA peak.")

main_args  = parser.add_argument_group('Main')
main_args.add_argument('-i','--input', dest='input', type=str, nargs='+',required=True,
                     help="A bam file, bam files, or single directory of bam files.")                    
main_args.add_argument('-o','--output',dest='out_root', type=str,required=True,
                     help="Name root to use for output files.")


config_args  = parser.add_argument_group('Config', 'Alter threhsolds, sizes, info tags e.t.c')
config_args.add_argument('--depth_threshold', dest='depth_threshold', type=int, default=10,
                    help="Need at least this many peaks to call an APA site. UNUSED")
config_args.add_argument('--region_size', dest='region_size', type=int, default=250,
                    help="Size upstream of APA to count reads. UNUSED")
config_args.add_argument('--minpolyA', dest='minpolyA', type=int, default=5,
                    help="Number of A to consider polyA. UNUSED")
config_args.add_argument('--cell_barcode_tag', dest='corrected_cell_barcode_tag', type=str, default='CB',
                     help="Corrected (exact-match) cell barcode bam tag used in bam input. UNUSED")
config_args.add_argument('--umi_tag', dest='corrected_umi_tag', type=str, default='UB',
                     help="Corrected (exact-match) UMI / molecular barcode bam tag used in bam input. UNUSED")
config_args.add_argument('--gene_tag', dest='gene_tag', type=str, default='GN',
                     help="Assigned gene barcode bam tag used in bam input. UNUSED")
                     


running_args  = parser.add_argument_group('Running', 'For changing how this script runs. Stop/start on polyA e.t.c')
running_args.add_argument('-p', '--peaks_gff', dest='peaks_gff', type=str,
                    help="If provided, use this gff file of peaks instead of making one from polyA reads. [DESCRIBE FORMAT]. UNUSED")
running_args.add_argument('-t', '--threads', dest='threads', type=int, default=1,
                    help="Num threads for multithreaded steps. UNUSED")
running_args.add_argument('--no_count', dest='skip_count', action='store_true', default=False,
                    help="Do not count reads in peaks. UNUSED" )                    
running_args.add_argument('--polyA_bams', dest='polyA_bams', nargs='*',
                    help="Skip polyA filtering step, and just use this already-processed polyA bam file(s) or directory of bams. UNUSED")


args = parser.parse_args()
print("Opitons: "+args+"\n")


###############################################################################
# MAIN
###############################################################################

def main (): 

    ## Fail early. 
    # Check tools in paths.
    # Check output clear
    
    
    ## Get and merge polyA bams
    input_bams = read_files_list_or_dir (args.input, filesuffix=".bam") 
    print("Processing input bam files: \n"+ "\n".join(input_bams) )
    # Check each bam.
    
    
    ## Get polyA peaks
    
    ## Count in polyA peaks
    
    ## Prep R object?


###############################################################################
# FUNCTIONS
###############################################################################







###############################################################################
# UTIL FUNCTIONS
###############################################################################


def read_files_list_or_dir (filesin, filesuffix=".bam") :
    
    actual_files = list()
    
    # If its a dir, just grab .bams and silently ignore .bai e.t.c
    if len(filesin) == 1 and os.path.isdir(filesin[0]) :
        for afile in os.listdir(filesin[0]) :
            if afile.endswith(filesuffix) and not os.path.isdir(afile) : 
                actual_files.append(os.path.join(filesin[0], afile))
    else :
        actual_files = filesin

    # Trust no one
    for afile in actual_files :
        if not os.path.exists(afile)       : sys.exit("File "+afile+" does not exist")                              
        if os.path.isdir(afile)            : sys.exit("Found directory "+afile+" in "+filesuffix+" file list (specify one dir OR list of files)")  
        if not afile.endswith(filesuffix)  : sys.exit("File "+afile+" does is not of expected type "+filesuffix)    
        if os.stat(afile).st_size == 0     : sys.exit("File "+afile+" is empty")                                   
    
    return actual_files


#def quick_bam_check (bam_file, cell_tag, umi_tag, gene_tag=None) :
    # Read chr names
    # Look for one of each of those tags in top n reads. 
    



###############################################################################
# GO...
###############################################################################

main()
