#!/usr/bin/env python3
import argparse
import os
import sys
import pysam
import subprocess
from collections import defaultdict

# ./polyApipe.py -i test_files/bams/mini1k.bam -o mini1k   #NB Small, and has no fwd polyA
# ./polyApipe.py -i test_files/bams_polyA/mini_polyA.bam -o xxxx   # alreayd polyA, but lots of r/f to see.


parser = argparse.ArgumentParser(description="Count reads in polyA peaks. "+
     "Given a set of annotated bams files (cell barcode, UMI code and gene) "+
     "(1) Find where peaks of polyA reads occur. "+
     "(2) Count how many reads per sample fall in each polyA peak.")

main_args  = parser.add_argument_group('Main')
main_args.add_argument('-i','--input', dest='input', type=str, nargs='+',required=True,
                     help="A bam file, bam files, or single directory of bam files.")                    
main_args.add_argument('-o','--output',dest='out_root', type=str,required=True,
                     help="Name root to use for output files.")

config_args  = parser.add_argument_group('PolyA peak thresholds', 'Peak-level configs')  
config_args.add_argument('--depth_threshold', dest='depth_threshold', type=int, default=10,
                    help="Need at least this many reads in a peak to call an APA site. UNUSED")
config_args.add_argument('--region_size', dest='region_size', type=int, default=250,
                    help="Size upstream of APA site to consider polyA reads. UNUSED")


config_args  = parser.add_argument_group('PolyA read thresholds', 'Read-level polyA filters e.t.c')
config_args.add_argument('--minMAPQ', dest='minMAPQ', type=int, default=10,
                    help="Minimum MAPQ mapping quality to consider.")      
config_args.add_argument('--minpolyA', dest='minpolyA', type=int, default=5,
                    help="Number of A to consider polyA.")
config_args.add_argument('--non_A_allowed', dest='nonA_allowed', type=int, default=0,
                    help="Number of non A bases permitted in polyA region (while still having --minpolyA As)")
config_args.add_argument('--misprime_A_count', dest='minpolyA', type=int, default=8,
                    help="Number of As seen in the last --misprime_in of aligned reads to label potential mispriming") 
config_args.add_argument('--misprime_in', dest='minpolyA', type=int, default=10,
                    help="Look for --misprime_A_count As in this many nucleotides at the end of a reaad alignment when labelling potential misprime")                     

                    
config_args  = parser.add_argument_group('Bam tags', 'For specifiing umi, cell, genes e.t.c')                                                      
config_args.add_argument('--cell_barcode_tag', dest='corrected_cell_barcode_tag', type=str, default='CB',
                     help="Corrected (exact-match) cell barcode bam tag used in bam input. UNUSED")
config_args.add_argument('--umi_tag', dest='corrected_umi_tag', type=str, default='UB',
                     help="Corrected (exact-match) UMI / molecular barcode bam tag used in bam input. UNUSED")
config_args.add_argument('--gene_tag', dest='gene_tag', type=str, default='GN',
                     help="Assigned gene barcode bam tag used in bam input. UNUSED")

                     


running_args  = parser.add_argument_group('Running', 'For changing how this script runs. Stop/start on polyA step e.t.c')
running_args.add_argument('-p', '--peaks_gff', dest='peaks_gff', type=str,
                    help="If provided, use this gff file of peaks instead of making one from polyA reads. [DESCRIBE FORMAT]. UNUSED")
running_args.add_argument('-t', '--threads', dest='threads', type=int, default=1,
                    help="Num threads for multithreaded steps. UNUSED")
running_args.add_argument('--no_count', dest='skip_count', action='store_true', default=False,
                    help="Do not count reads in peaks. UNUSED" )                    
running_args.add_argument('--polyA_bams', dest='polyA_bams', nargs='*',
                    help="Skip polyA filtering step, and just use this already-processed polyA bam file(s) or directory of bams. UNUSED")


args = parser.parse_args()

#print("Opitons: "+args+"\n")


###############################################################################
# MAIN
###############################################################################


        
        


def main (): 

    ## Fail early. 
    # parameters
    check_params_ok(args)
    # Check tools in paths.
    check_tools_available()
    # Check output clear
    
    polyA_bam_root = args.out_root+"_polyA"


    
    ## Get and merge polyA bams
    input_bams = read_files_list_or_dir (args.input, filesuffix=".bam") 
    print("Finding input bam files: \n"+ "\n".join(input_bams) )
    print("")
    
    # Check each bam is ok
    print("Checking each bam:")
    for input_bam in input_bams :
        quick_bam_check (input_bam, args.corrected_cell_barcode_tag, args.corrected_umi_tag, args.minMAPQ)
    print("")
    
    ## Get polyA reads
    print("Finding polyA reads in each input bam file: ")
    
    
    
    process_bams_to_polyA_bam(input_bams, polyA_bam_root, args.minpolyA, args.minMAPQ, args.nonA_allowed) 
    print("")
        
    ## Count in polyA peaks




###############################################################################
# FUNCTIONS - polyA bams
###############################################################################

def process_bams_to_polyA_bam (input_bams, polyA_bam_root, minpolyA, minMAPQ, nonA_allowed)  :

    polyA_bam_file =  polyA_bam_root+".bam"
    polyA_bam_dir  =  polyA_bam_root+"_individual_bams"

    if len(input_bams) == 1 :
        make_only_polyA_bam(input_bams[0], polyA_bam_file, minpolyA, minMAPQ, nonA_allowed)
    
    else :  #If more than one, then dump em in separately and merge.
        
        try:  
            os.mkdir(polyA_bam_dir)
        except OSError:  
            #sys.exit("Couldn't create polyA output directory "+polyA_bam_dir)
            print("PolyA output dir exists arleady (TEMP)")
            pass
        
        #Potentially pararllelise here?
    
        # Filter down  each bam file to just polyA reads
        polyA_bamfile_inds = list()
        for input_bam in input_bams :
            polyA_bamfile_ind = os.path.join(polyA_bam_dir, os.path.basename(input_bam))
            polyA_bamfile_inds.append(polyA_bamfile_ind )
            make_only_polyA_bam(input_bam, polyA_bamfile_ind, minpolyA, minMAPQ, nonA_allowed)
            pysam.index(polyA_bamfile_ind)  # Samtools index (was indexed already :. arleady sorted.)
            print("Got polyA reads from"+input_bam)
            
        # Now merge result.
        print(polyA_bam_file+" ".join(polyA_bamfile_inds))
        merge_params = [polyA_bam_file] + polyA_bamfile_inds
        pysam.merge( *merge_params )

    pysam.index(polyA_bam_file)
    if not os.path.exists(polyA_bam_file+".bai")       : sys.exit("No bam index made for polyA-only bam file "+polyA_bam_file+" Something broke.")       
    print("Got all polyA reads into "+polyA_bam_file)





def make_only_polyA_bam (bam_file, polyA_bamfile, minpolyA, minMAPQ, nonA_allowed) :
    
    
    bam      = pysam.AlignmentFile(bam_file,'rb')
    polyAbam = pysam.AlignmentFile(polyA_bamfile, "wb", template=bam)

    for read in bam.fetch(until_eof=True) :
        
        strand  = "-" if read.is_reverse else "+"
        readseq = read.query_sequence     
        #print("strand=\t"+strand)
        #print("mapq=\t"+str(read.mapping_quality))
        #print("readseq=\t"+readseq)   
        #print("cigar=\t"+read.cigarstring)
        #print(str(read.query_alignment_start)+" - "+str(read.query_alignment_end) + "("+strand+") == "+str(read.query_length)+"")
        
                
        if read.mapping_quality < minMAPQ : continue

        is_polyA = False
        if strand == "+" :
            
            soft_match_size = read.query_length - read.query_alignment_end - 1 # 0-based                     
            if soft_match_size >= minpolyA :
                soft_match_seq = read.query_sequence[-soft_match_size:]
                num_As         = soft_match_seq.upper().count('A')
                
                if num_As >= minpolyA  and  num_As >= ( len(soft_match_seq) - nonA_allowed ) :
                    is_polyA = True
                    
            
        else : # - strand reads have TTTTTT at the alignment start
            soft_match_size = read.query_alignment_start  # 0-based           
            if soft_match_size >= minpolyA :
                soft_match_seq = read.query_sequence[0:soft_match_size]
                num_As         = soft_match_seq.upper().count('T')

                if num_As >= minpolyA  and  num_As >= ( len(soft_match_seq) - nonA_allowed ) :
                    is_polyA = True

        
        if is_polyA : 
            polyAbam.write(read) 
    
    polyAbam.close()    
    bam.close()


###############################################################################
# FUNCTIONS - polyA peaks
###############################################################################

#def process_polyA_ends_to_peaks(polyA_bams, depth_threshold, region_size) : 








###############################################################################
# UTIL FUNCTIONS
###############################################################################

def check_params_ok (args) :
    if args.minpolyA < 0 : 
        sys.exit("--minpolyA should be a positive integer")
    if args.region_size < 0 : 
        sys.exit("--region_size should be a positive integer")
    if args.depth_threshold < 0 : 
        sys.exit("--depth_threshold should be a positive integer")    
    if args.nonA_allowed < 0 : 
        sys.exit("--non_A_allowed should be a positive integer")
    

        
def check_outputs_clean () :
    
    # out root not / . remove trainiling /
    pass




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
    seenit = dict()
    for afile in actual_files :
        if not os.path.exists(afile)         : sys.exit("File "+afile+" does not exist")                              
        if os.path.isdir(afile)              : sys.exit("Found directory "+afile+" in "+filesuffix+" file list (specify one dir OR list of files)")  
        if not afile.endswith(filesuffix)    : sys.exit("File "+afile+" does is not of expected type "+filesuffix)    
        if os.stat(afile).st_size == 0       : sys.exit("File "+afile+" is empty")                                   
        if os.path.basename(afile) in seenit : sys.exit("Filename '"+afile+"' is duplicated in input (all input bamfiles must have unique names)")
        seenit[os.path.basename(afile)] = 1 # now seenit
        
    return actual_files



def check_tools_available () :
    print("Checking for tools:")
    try:
        subprocess.call(["samtools","--version"])
    except OSError as e:
        sys.exit("Could not find samtools in PATH. The samtools package should be installed and in PATH.")
    print("Tools ok\n")



def quick_bam_check (bam_file, cell_tag, umi_tag, minMAPQ) :
    # is there an index?
    bam_index = bam_file+".bai"
    
    if not os.path.exists(bam_index)       : sys.exit("No bame index "+bam_index+" for file "+bam+" Use samtools index input bams.")       

    # Read chr names  (actually, only from bam, so no need.)
    # Look for one of each of those tags in top n reads. Counting won't happen without them!
    top_n = 1000
    n     = 0
    seen_cell_tag = False
    seen_umi_tag  = False
    seen_map_qual_ok = False
    #seen_gene_tag = False Actually, don't check. Quitely likely there's no GN near the start, and it works without anywa.y.
    bam = pysam.AlignmentFile(bam_file,'rb')

    
    for read in bam.fetch(until_eof=True) :
        n = n+1 
        if (n > top_n) : 
            break
            
        if read.mapping_quality >= minMAPQ : seen_map_qual_ok = True
        
        try:
            cb = read.get_tag(cell_tag)
            seen_cell_tag = True
        except KeyError:
            pass
        try:
            cb = read.get_tag(umi_tag)
            seen_umi_tag = True
        except KeyError:
            pass

    bam.close()

    if ( not seen_cell_tag or not seen_umi_tag ) :
        
        missing = cell_tag if not seen_cell_tag else "" 
        missing = (missing +" "+ umi_tag) if not seen_umi_tag else missing 
        
        sys.exit( "Checked the first "+str(top_n)+" reads of "+bam_file+" and did not see any "+missing+ " tags.\n" +
        "Fix by annotating cell barcodes and/or UMIs in tags in bam file, and specifying the tags with"+
        " (--cell_barcode_tag / --umi_tag)") 

    if (not seen_map_qual_ok) :
        sys.exit( "Checked the first "+str(top_n)+" reads of "+bam_file+" and did not see any reads that have min map quality "+str(minmapQ))


    print (bam_file+" ok")

        
    



###############################################################################
# GO...
###############################################################################

main()
