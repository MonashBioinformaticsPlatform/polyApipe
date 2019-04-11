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
                    help="Need at least this many reads in a peak to call an APA site.")
config_args.add_argument('--region_size', dest='region_size', type=int, default=250,
                    help="Size upstream of APA site to consider polyA reads.")


config_args  = parser.add_argument_group('PolyA read thresholds', 'Read-level polyA filters e.t.c')
config_args.add_argument('--minMAPQ',          dest='minMAPQ',          type=int, default=10,
                    help="Minimum MAPQ mapping quality to consider.")      
config_args.add_argument('--minpolyA',         dest='minpolyA',         type=int, default=5,
                    help="Number of A to consider polyA.")
config_args.add_argument('--non_A_allowed',    dest='nonA_allowed',     type=int, default=0,
                    help="Number of non A bases permitted in polyA region (while still having --minpolyA As)")
config_args.add_argument('--misprime_A_count', dest='misprime_A_count', type=int, default=8,
                    help="Number of As seen in the last --misprime_in of aligned reads to label potential mispriming") 
config_args.add_argument('--misprime_in',      dest='misprime_in',      type=int, default=10,
                    help="Look for --misprime_A_count As in this many nucleotides at the end of a reaad alignment when labelling potential misprime")                     

                    
config_args  = parser.add_argument_group('Bam tags', 'For specifiing umi, cell, genes e.t.c')                                                      
config_args.add_argument('--cell_barcode_tag', dest='corrected_cell_barcode_tag', type=str, default='CB',
                     help="Corrected (exact-match) cell barcode bam tag used in bam input.")
config_args.add_argument('--umi_tag', dest='corrected_umi_tag', type=str, default='UB',
                     help="Corrected (exact-match) UMI / molecular barcode bam tag used in bam input.")
config_args.add_argument('--gene_tag', dest='gene_tag', type=str, default='GN',
                     help="Assigned gene barcode bam tag used in bam input.")

                     


running_args  = parser.add_argument_group('Running', 'For changing how this script runs. Stop/start on polyA step e.t.c')
running_args.add_argument('--no_peaks', dest='skip_peaks', action='store_true', default=False,
                    help="Stop after making polyA bams. Do not try to find peaks in polyA files (implies --no_count)" )            
running_args.add_argument('--no_count', dest='skip_count', action='store_true', default=False,
                    help="Do not count reads in peaks." )                    
running_args.add_argument('--polyA_bams', dest='polyA_bams', action='store_true', default=False,
                    help="Skip polyA filtering step, the bams specified with '-i' are already filtered to polyA-containing reads only.")
running_args.add_argument('-p', '--peaks_gff', dest='peaks_gff', type=str, default=None,
                    help="If provided, use this gff file of peaks instead of making one from polyA reads. Will still try to make those polyA bams unless --polyA_bams also specified. [DESCRIBE FORMAT]. ")
running_args.add_argument('-t', '--threads', dest='threads', type=int, default=1,
                    help="Num threads for multithreaded steps. UNUSED")

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
    
    polyA_bam_root   = args.out_root+"_polyA"
    polyA_peaks_gff  = args.out_root+"_polyA_peaks.gff"
    polyA_bam        = polyA_bam_root+".bam"
            
    ## Get and merge polyA bams
    input_bams = read_files_list_or_dir (args.input, filesuffix=".bam") 
    print("Finding input bam files: \n"+ "\n".join(input_bams) )

    
    # Check each bam is ok
    print("\nChecking each bam:")
    for input_bam in input_bams :
        quick_bam_check (input_bam, args.corrected_cell_barcode_tag, args.corrected_umi_tag, args.minMAPQ)

    
    ## Get polyA reads
    if not args.polyA_bams :   
        print("\nFinding polyA reads in each input bam file: ")
        process_bams_to_polyA_bam(  input_bams     = input_bams, 
                                    polyA_bam_root = polyA_bam_root, 
                                    minpolyA       = args.minpolyA, 
                                    minMAPQ        = args.minMAPQ, 
                                    nonA_allowed   = args.nonA_allowed) 
                                    
    else : # Or heres one we prepared earlier
        print("\nThe bam files provided have already been filtered to polyA reads. Skipping polyA filter step.")   
        if len(input_bams) == 1  :
            polyA_bam = input_bams[0]
        else :
            print("\nMerging polyA bam files.")
            merge_bams(input_bams, polyA_bam, indexit=True)
                

    if args.skip_peaks : sys.exit( "\nRequested no peak calling. Finished." )                
                                
        
    ## Define polyA reads
    if not args.peaks_gff :
        print("\nFinding peaks in merged polyA file: ")
        process_polyA_ends_to_peaks(polyA_bam         = polyA_bam , 
                                    polyA_peaks_gff   = polyA_peaks_gff,
                                    depth_threshold   = args.depth_threshold, 
                                    region_size       = args.region_size,
                                    corrected_cell_barcode_tag = args.corrected_cell_barcode_tag, 
                                    corrected_umi_tag = args.corrected_umi_tag, 
                                    gene_tag          = args.gene_tag,
                                    misprime_A_count  = args.misprime_A_count,
                                    misprime_in       = args.misprime_in)
    else :
       polyA_peaks_gff = args.peaks_gff
       print("\n Using provided gff file: " + polyA_peaks_gff )
       
        
 
    if args.skip_count : sys.exit( "\nRequested no counting. Finished." )       
 
 
 
 
    ## Count in polyA peaks
    print("\nCount reads in polyA peaks: ")
    

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
        merge_bams (polyA_bamfile_inds, polyA_bam_file)

    pysam.index(polyA_bam_file)
    if not os.path.exists(polyA_bam_file+".bai")       : sys.exit("No bam index made for polyA-only bam file "+polyA_bam_file+" Something broke.")       
    print("Got all polyA reads into "+polyA_bam_file)





def make_only_polyA_bam (bam_file, polyA_bamfile, minpolyA, minMAPQ, nonA_allowed) :
    # SLIGHTLY DIFFERETENT TO only-polA.pl!!!
    # option to allow n nonA chars in the soft maasked portion (assuming suffient minpolyA As.)
    # Rather than allowing nonA chars at the very end (ie stray adaptor) soft masked portion (assuming suffient minpolyA As.)
    # May / may not matter...
    
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

def process_polyA_ends_to_peaks(polyA_bam, polyA_peaks_gff, depth_threshold, region_size,
    corrected_cell_barcode_tag, corrected_umi_tag, gene_tag,
    misprime_A_count, misprime_in) : 

    # Some rough numbers to flag where things go wrong.
    reads_total       = 0 
    reads_considered  = 0
    reads_having_gene = 0
    
    
    samfile = pysam.AlignmentFile(polyA_bam,'rb')
    #ends[contig][pos][strand]
    ends = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    misprime_site = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    genes = defaultdict(set)


    print("Processing read ends", file=sys.stderr)
    for read in samfile.fetch(until_eof=True):
        reads_total +=1
        #print(read)
        try:
          cb = read.get_tag(corrected_cell_barcode_tag)
          ub = read.get_tag(corrected_umi_tag)
          reads_considered += 1
        except KeyError:
          # Skip reads without cell or umi
          continue

        strand = "-" if read.is_reverse else "+"

        # Try to avoid "mis-primed" reads. 
        # If the aligned end of read has a number of A's, tag as potentionally mis-primed
        seq = read.query_alignment_sequence
        if read.is_reverse:
          numA = seq[:misprime_in].count('T')
        else:
          numA = seq[-misprime_in:].count('A')
        misprime = numA>=misprime_A_count

        pos = read.reference_start if read.is_reverse else (read.reference_end-1)

        try:
          gn = read.get_tag(gene_tag)
          reads_having_gene +=1
        except KeyError:
          gn=None

        # Ignore reads tagged as duplicates
        if read.is_duplicate:
          continue

        if gn:   # POTENTIALLY PROBLEMATIC FORMATTING CUSTOMISATION
          # e.g. RHOT2 or RHOF;LINC01089
          genes[read.reference_name+str(pos)+strand].update( gn.split(";") )

        # Count depth for end
        ends[read.reference_name][pos][strand] += 1

        if misprime:
          misprime_site[read.reference_name][pos][strand] += 1

    # regions[contig][start] = (end,strand)
    regions = defaultdict(lambda: defaultdict(int))

    # Create potential counting regions
    #ends[contig][pos][strand]
    print("Creating regions", file=sys.stderr)
    for contig in sorted(ends.keys()):
      for pos in sorted(ends[contig].keys()):
        for strand in sorted(ends[contig][pos].keys()):
          count = ends[contig][pos][strand]
          misprime = misprime_site[contig][pos][strand] >= count/2         # If at least half the reads are tagged as misprimed (I think it should be all, or none)
          if count>=depth_threshold:
            names=genes[contig+str(pos)+strand]   
            if strand=='+':
              regions[contig][max(pos-region_size,0)] = (count,pos,strand,names,misprime)
            else:
              regions[contig][pos] = (count,pos+region_size,strand,names,misprime)

    # Find regions "too close" and drop the one with lowest depth
    print("Handling overlapping regions", file=sys.stderr)
    for contig in sorted(regions.keys()):
      kept = {}
      for strand in ['+','-']:
        starts = sorted(x for x in regions[contig].keys() if regions[contig][x][2]==strand)
        for i in range(len(starts)):
          start = starts[i]
          (depth, end, _, names, misprime) = regions[contig][start]

          best = True
          # Check if there was a better region upstream
          i2 = i-1
          while (i2>=0 and best):
            (depth2, end2, _, _, misprime2) = regions[contig][starts[i2]]
            if end2<start:
              break              # Doesn't overlap

            if misprime==misprime2 and depth2>depth:
              best = False       # Other one is better (deeper)
            elif misprime and not misprime2:
              best = False       # Other one is not misprimed

            i2 -= 1

          # Check if there is a better region downstream
          i2=i+1
          while (i2<len(starts) and best):
            (depth2, end2, _, _, misprime2) = regions[contig][starts[i2]]
            if end<starts[i2]:
              break              # Doesn't overlap

            if misprime==misprime2 and depth2>=depth:
              best = False       # Other one is better (deeper)
            elif misprime and not misprime2:
              best = False       # Other one is not misprimed

            i2 += 1

          # We're the best of the overlapping!
          if best:
            kept[start] = (depth, end, strand, names, misprime)
      regions[contig] = kept

    # Print regions
    print("Writing output "+polyA_peaks_gff, file=sys.stderr)
    gff = open(polyA_peaks_gff, "w")
    gff.write("##gff\n") ##gff-version 3

    for contig in sorted(regions.keys()):
      for pos in sorted(regions[contig].keys()):
        (depth,end,strand,names,misprime) = regions[contig][pos]
    
    
        # 0 to 1 based (was bed?)
        pos = pos + 1 if strand == "-" else pos + 2
        end = end + 1 if strand == "+" else end  

    
        gene_name=""
        if len(names)==0:
          gene_name = "Unknown"  
        else : 
          gene_name = ",".join(names)
        str_for_name = "f" if strand == "+" else "r"
        pos_for_name = pos if strand == "-" else end
        peak_name    = gene_name+":%s_%d_%s"%(contig, pos_for_name, str_for_name)
    

        # There is not bed, only gff
        #12	polyAends	polyApeak	121802524	121802774	.	-	.	peakgene="RHOF"; peak="RHOF:121779162"; peakdepth="20"; misprime="False";    
        detail='peakgene="%s"; peak="%s"; peakdepth="%d"; misprime="%s";'%(gene_name, peak_name, depth, str(misprime) )
        gff.write("\t".join(str(x) for x in [contig,"polyAends","polyAends",pos,end,".",strand,".", detail,"\n"  ]))
        
    
    gff.close()
    
    # Summary
    print("Processed %d polyA reads, %d included (had %s and %s tags), %d of which had a tagged gene (0 is ok if no genes anotated with %s)"%(
          reads_total, reads_considered,corrected_cell_barcode_tag ,corrected_umi_tag, reads_having_gene, gene_tag) )
    
    


###############################################################################
# FUNCTIONS - counting
###############################################################################



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

        
    
def merge_bams (bams_in, merged_bam, indexit=False):
    merge_params = [merged_bam] + bams_in
    pysam.merge( *merge_params )
    if indexit :
        pysam.index(merged_bam)
    


###############################################################################
# GO...
###############################################################################

main()
