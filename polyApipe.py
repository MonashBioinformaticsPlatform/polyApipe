#!/usr/bin/env python3
import argparse
import os
import sys
import pysam
import subprocess
import warnings
from collections import defaultdict
from distutils.version import StrictVersion

# module load subread

# ./polyApipe.py -i test_files/demo/  -o xxxx
# ./polyApipe.py -i test_files/demo/SRR5259354_demo.bam  -o xxxxSRR5259354_demo
# ./polyApipe.py -i test_files/demo/SRR5259354_demo.bam test_files/demo/SRR5259422_demo.bam  -o xxxxtwo
# ./polyApipe.py -i data/demo/ -o polyApiper/inst/extdata/demo_dataset/demo
# ./polyApipe.py -i data/demo/SRR5259422_demo.bam  -o xxSRR5259422_demo
# ./polyApipe.py -i test_files/bams_polyA/mini_polyA.bam -o xxxx   # alreayd polyA, but lots of r/f to see.


__version__='0.1.0' # Keep in sync with polyapiper

parser = argparse.ArgumentParser(
    description="Count reads in polyA peaks. "+
        "Given a set of bams files (with cell barcode and UMI code tags) "+
        "(1) Find where peaks of polyA reads occur. "+
        "(2) Count how many reads per sample fall in each polyA peak.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter )#

parser.add_argument('-v', '--version', action='version',
                    version='%(prog)s {version}'.format(version=__version__))


main_args  = parser.add_argument_group('Main')
main_args.add_argument('-i','--input', dest='input', type=str, nargs='+',required=True,
                     help="A bam file, bam files, or single directory of bam files.")                    
main_args.add_argument('-o','--output',dest='out_root', type=str,required=True,
	             help="Name root to use for output files.")

config_args  = parser.add_argument_group('PolyA peak options', 'Set peak-level configs')
config_args.add_argument('--depth_threshold', dest='depth_threshold', type=int, default=10,
                    help="Need at least this many reads in a peak to call an APA site.")
config_args.add_argument('--region_size', dest='region_size', type=int, default=250,
                    help="Size upstream of APA site to consider polyA reads.")


config_args  = parser.add_argument_group('PolyA read options', 'Set read-level configs for defining polyA reads')
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

                    
config_args  = parser.add_argument_group('Bam tags', 'For specifiing umi, cell e.t.c')
config_args.add_argument('--cell_barcode_tag', dest='corrected_cell_barcode_tag', type=str, default='CB',
                     help="Corrected (exact-match) cell barcode bam tag used in bam input.")
config_args.add_argument('--umi_tag', dest='umi_tag', type=str, default='UR',
                     help="Uncorrected UMI / molecular barcode bam tag used in bam input. May contain mismatches.")               


running_args  = parser.add_argument_group('Running', 'For changing how this script runs. Stop/start on polyA step e.t.c')
#running_args.add_argument('-v', '--version', dest='print_version', action='store_true', default='False',
#                    help="Print version and exit.")
running_args.add_argument('-t', '--threads', dest='threads', type=int, default=1,
                    help="Num threads for multithreaded steps.")
running_args.add_argument('--no_peaks', dest='skip_peaks', action='store_true', default=False,
                    help="Stop after making polyA bams. Do not try to find peaks in polyA files (implies --no_anno --no_count)" )          
running_args.add_argument('--no_count', dest='skip_count', action='store_true', default=False,
                    help="Stop after making merged polyA peaks gff file. Do not try to count them, or annotate bams with them." )
running_args.add_argument('--polyA_bams', dest='polyA_bams', action='store_true', default=False,
                    help="Skip polyA filtering step, the bams specified with '-i' are already filtered to polyA-containing reads only.")            
running_args.add_argument('--peak_anno_bams', dest='peak_anno_bams', action='store_true', default=False,
                    help="DEPRECATED - The bams provided have already been labelled with peaks regions e.t.c Jump to immediate counting.")
running_args.add_argument('--peaks_gff', dest='peaks_gff', type=str, default=None,
                    help="If provided, use this gff file of peaks instead of making one from polyA reads. "+
                         "Skips polyA filtering steps, incompatable with --polyA_bams "+
                         "This is a gtf format specifically as output by this script. See example data.")
running_args.add_argument('--keep_interim_files', dest='keep_interim_files', action='store_true', default=False,
                    help="Don't delete the intermediate files (merged polyA, annotated input bams)." +
                         "For debugging or piecemeal runs." )
running_args.add_argument('--skipchecks', dest='skip_checks', action='store_true', default=False,
                    help=argparse.SUPPRESS) # nothing to see here, move along.


args = parser.parse_args()

###############################################################################
# MAIN
###############################################################################


        
        


def main (): 

    ## Fail early. 
    # parameters
    check_params_ok(args)
    # Check tools in paths.
    if not args.skip_checks : check_tools_available()

    # Check output clear
    
    polyA_bam_root     = args.out_root+"_polyA"
    polyA_bam          = polyA_bam_root+".bam"
    polyA_peaks_gff    = args.out_root+"_polyA_peaks.gff"
    peak_anno_bam_root = args.out_root+"_peakanno"
    counts_root        = args.out_root+"_counts"
    
    
    
            
    ## Get and merge polyA bams
    input_bams     = read_files_list_or_dir (args.input, filesuffix=".bam") 
    input_from_dir = os.path.isdir(args.input[0]) or len(input_bams) > 1 # Already sanity checked, if dir, only one dir
    print("Finding input bam files: \n"+ "\n".join(input_bams) )


    
    # Check each bam is ok
    print("\nChecking each bam:")
    for input_bam in input_bams :
        quick_bam_check (input_bam, args.corrected_cell_barcode_tag, args.umi_tag, args.minMAPQ)


    ## Define polyA peaks
    if not args.peaks_gff :
    
        ## Get polyA reads
        if not args.polyA_bams :
            print("\nFinding polyA reads in each input bam file: ")
            process_bams_to_polyA_bam(  input_bams         = input_bams,
                                        polyA_bam_root     = polyA_bam_root,
                                        minpolyA           = args.minpolyA,
                                        minMAPQ            = args.minMAPQ,
                                        nonA_allowed       = args.nonA_allowed)

        else : # Or heres one we prepared earlier
            print("\nThe bam files provided have already been filtered to polyA reads. Skipping polyA filter step.")
            if len(input_bams) == 1  :
                polyA_bam = input_bams[0]
            else :
                print("\nMerging polyA bam files.")
                merge_bams(input_bams, polyA_bam, indexit=True)


        if args.skip_peaks :
            print( "\nRequested no peak calling. Finished." )
            sys.exit(0)

        

        print("\nFinding peaks in merged polyA file: ")
        process_polyA_ends_to_peaks(polyA_bam          = polyA_bam ,
                                    polyA_peaks_gff    = polyA_peaks_gff,
                                    depth_threshold    = args.depth_threshold,
                                    region_size        = args.region_size,
                                    corrected_cell_barcode_tag = args.corrected_cell_barcode_tag, 
                                    umi_tag            = args.umi_tag,
                                    misprime_A_count   = args.misprime_A_count,
                                    misprime_in        = args.misprime_in)

        # Done peak finding, don't need merged polyA (if any) bam anymore (unless polyA bam(s) supplied)
        if not args.keep_interim_files and not args.polyA_bams:
            if input_from_dir:
                os.remove(polyA_bam)
                os.remove(polyA_bam + ".bai")

    else :
        # ALready had peaks file, so no need for polyA stuff at all.
        polyA_peaks_gff = args.peaks_gff
        print("\n Using provided gff file: " + polyA_peaks_gff )


    if args.skip_count :
        print( "\nRequested no counting. Finished.")
        sys.exit(0)
 

 
    ## Annotate ORIGINAL bams with peaks
    print("\nMaking filtered bam files annotated with peak hits (XT tag):")
    anno_bams = None
    if (not args.peak_anno_bams) :
        anno_bams = make_peak_hits_annotated_bams ( 
                                    input_bams          = input_bams,
                                    input_from_dir      = input_from_dir,
                                    polyA_peaks_gff     = polyA_peaks_gff,
                                    peak_anno_bam_root  = peak_anno_bam_root,
                                    corrected_cell_barcode_tag = args.corrected_cell_barcode_tag,
                                    threads             = args.threads)

    else :
        anno_bams = input_bams
    

    ## Now count with umi tools.
    print("\nCounting reads per peak per cell:")
    count_from_annotated_bams(      peak_anno_bams = anno_bams, 
                                    input_from_dir = input_from_dir, 
                                    counts_root    = counts_root, 
                                    corrected_cell_barcode_tag = args.corrected_cell_barcode_tag, 
                                    umi_tag        = args.umi_tag )

    # Remove those large anno bam files (unless they were the input!)
    if not args.keep_interim_files and not args.peak_anno_bams :
        for anno_bam in anno_bams :
            os.remove(anno_bam)
            os.remove(anno_bam+".bai")
        if input_from_dir :
            os.rmdir(peak_anno_bam_root)

 
    print("\nDone!")

    
    
    
    
    

###############################################################################
# FUNCTIONS - polyA bams
###############################################################################

def process_bams_to_polyA_bam (input_bams, polyA_bam_root, minpolyA, minMAPQ, nonA_allowed)  :

    polyA_bam_file     =  polyA_bam_root+".bam"
    polyA_bam_dir      =  polyA_bam_root+"_individual_bams"


    # Either the input was a single file, or a directory with a single file. 
    # Either way, no need for an individual_bams directory.
    if len(input_bams) == 1 :
        make_only_polyA_bam(input_bams[0], polyA_bam_file, minpolyA, minMAPQ, nonA_allowed)
    
    else :  #If more than one, then dump em in separately and merge.
        
        try:  
            os.mkdir(polyA_bam_dir)
        except OSError:  
            #sys.exit("Couldn't create polyA output directory (already exists?) "+polyA_bam_dir)
            print("PolyA output dir exists already (TEMP)")
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
    # SLIGHTLY DIFFERETENT TO only-polyA.pl!!!
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
    corrected_cell_barcode_tag, umi_tag, misprime_A_count, misprime_in) :

    # Some rough numbers to flag where things go wrong.
    reads_total       = 0 
    reads_considered  = 0

    ends          = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    misprime_site = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    print("Processing read ends", file=sys.stderr)
    samfile = pysam.AlignmentFile(polyA_bam,'rb')

    for read in samfile.fetch(until_eof=True):
        reads_total +=1

        try:
          cb = read.get_tag(corrected_cell_barcode_tag)
          ub = read.get_tag(umi_tag)
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


        # Ignore reads tagged as duplicates
        if read.is_duplicate:
            continue

        # Count depth for end
        ends[read.reference_name][pos][strand] += 1

        if misprime:
            misprime_site[read.reference_name][pos][strand] += 1

    # regions[contig][start] = (count, pos, strand, misprime)
    regions = defaultdict(lambda: defaultdict(int))

    # Create potential counting regions
    # ends[contig][pos][strand]
    print("Creating regions", file=sys.stderr)
    for contig in sorted(ends.keys()):
        for pos in sorted(ends[contig].keys()):
            for strand in sorted(ends[contig][pos].keys()):
                count    = ends[contig][pos][strand]
                misprime = misprime_site[contig][pos][strand] >= count/2  # If at least half the reads are tagged as misprimed (I think it should be all, or none)
                if count>=depth_threshold:
                    if strand == '+':
                        regions[contig][max(pos-region_size,0)] = (count, pos, strand, misprime)
                    else:
                        regions[contig][pos] = (count, pos+region_size, strand, misprime)

    # Find regions "too close" and drop the one with lowest depth
    print("Handling overlapping regions", file=sys.stderr)
    for contig in sorted(regions.keys()):
        kept = {}
        for strand in ['+','-']:
            starts = sorted(x for x in regions[contig].keys() if regions[contig][x][2]==strand)
            for i in range(len(starts)):
                start = starts[i]
                (depth, end, _, misprime) = regions[contig][start]

                best = True
                # Check if there was a better region upstream
                i2 = i-1
                while (i2>=0 and best):
                    (depth2, end2, _, misprime2) = regions[contig][starts[i2]]
                    if end2<start:
                        break              # Doesn't overlap

                    if depth2>depth:
                        best = False       # Other one is better (deeper)
                    i2 -= 1

                # Check if there is a better region downstream
                i2=i+1
                while (i2<len(starts) and best):
                    (depth2, end2, _, misprime2) = regions[contig][starts[i2]]
                    if end<starts[i2]:
                        break              # Doesn't overlap

                    if depth2>=depth:
                        best = False       # Other one is better (deeper)
                    i2 += 1

                # We're the best of the overlapping!
                if best:
                    kept[start] = (depth, end, strand, misprime)
        regions[contig] = kept

    # Print regions
    print("Writing output "+polyA_peaks_gff, file=sys.stderr)
    gff = open(polyA_peaks_gff, "w")
    gff.write("##gff\n") ## atm its a gtf / gff2 like ensembl


    for contig in sorted(regions.keys()):
        for pos in sorted(regions[contig].keys()):

            (depth,end,strand,misprime) = regions[contig][pos]

            # 0 to 1 based (was bed?)
            pos = pos + 1 if strand == "-" else pos + 2
            end = end + 1 if strand == "+" else end
            str_for_name = "f" if strand == "+" else "r"
            pos_for_name = pos if strand == "-" else end
            peak_name    = "%s_%d_%s"%(contig, pos_for_name, str_for_name)

            # There is no bed, only gtf
            detail='peak="%s"; peakdepth="%d"; misprime="%s";'%(peak_name, depth, str(misprime) )
            gff.write("\t".join(str(x) for x in [contig,"polyAends","polyAends",pos,end,".",strand,".", detail ])+"\n")
        
    
    gff.close()



    # Summary
    print("Processed %d polyA reads, %d included (had %s and %s tags)"%(
          reads_total, reads_considered,corrected_cell_barcode_tag ,umi_tag) )
    



###############################################################################
# FUNCTIONS - counting
###############################################################################

def make_peak_hits_annotated_bams (input_bams, input_from_dir, polyA_peaks_gff, peak_anno_bam_root ,  corrected_cell_barcode_tag,  threads) :
    
    
     ## Figure out the output.
    peak_anno_bams = list()
    if not input_from_dir :
        peak_anno_bams.append(peak_anno_bam_root+".bam")
         
    else :  
        
        try:  
            os.mkdir(peak_anno_bam_root)
        except OSError:  
            #sys.exit("Couldn't create peak annotated bam directory (already exists?)"+peak_anno_bam_dir)
            print("Peak annotated bam dir exists already (TEMP)")
            pass
            
        for input_bam in input_bams : 
            peak_anno_bams.append(os.path.join(peak_anno_bam_root, os.path.basename(input_bam)))


    # Process each.
    for n in range(0,len(input_bams)) : 
        input_bam     = input_bams[n]
        peak_anno_bam = peak_anno_bams[n] 
        
        print(input_bam +" to "+ peak_anno_bam)
        make_peak_hit_annotated_bam (input_bam, polyA_peaks_gff , peak_anno_bam, corrected_cell_barcode_tag, threads)
    
    
    return(peak_anno_bams)



#uses UR?
def make_peak_hit_annotated_bam (input_bam, peaks_gff , peak_anno_bam, corrected_cell_barcode_tag, threads):
    
    #NB: multiple runs on the same bams with different annos could clobber each other if simultaneous!
    # Due to how featurecounts names bams. (-o is only the counts table, not bams)
    #  => temp feature counts directory based on feature name to avoid this.

    temp_feature_counts_dir    = os.path.basename(peak_anno_bam)+"_"+os.path.basename(peaks_gff)+"_temp"
    try:  
        os.mkdir(temp_feature_counts_dir)
    except OSError:  
        #sys.exit("Couldn't create temporary dir for featureCounts run (already exists?)"+temp_feature_counts_dir)
        print("Couldn't create temporary dir for featureCounts run (already exists?)"+temp_feature_counts_dir)


    # Run featureCounts
    # Adds XT tag to bams:  XT:Z:1_16442_r
    # -f count at local feature leve, not metafeature (gene)
    print(input_bam + " running featureCounts")
    temp_unfiltered            = os.path.join(temp_feature_counts_dir, os.path.splitext(os.path.basename(peak_anno_bam))[0]) # is an outputfile 
    temp_unfiltered_bam        = os.path.join(temp_feature_counts_dir, os.path.basename(input_bam )+ ".featureCounts.bam") # bam named via input file, not -o
    temp_unfiltered_summary    = temp_unfiltered + ".summary"
    #saved_unfiltered_summary   = os.path.basename(temp_unfiltered) + "_count_totals.txt" #+"."+os.path.basename(peak_anno_bam)
    
    fc_cmd  = ["featureCounts", 
               input_bam, 
               "-t", "polyAends", "-g", "peak", "-F", "GTF", "-f", "-s", "1", "-R", "BAM",
               "-T", str(threads),
               "-a", peaks_gff,
               "-o", temp_unfiltered  ]
    ranok = subprocess.call( fc_cmd ) 
    if not ranok == 0 :
        sys.exit("Failed to run featureCounts correctly with cmd\n"+ " ".join(fc_cmd))


    ## Filter out everythign without a CB tag (else umitools complains) - UR doesn't matter here.
    # Also, filter out annthing without an annotation, as this dramamtically reduces the size of the bam (which aren't being kept anyway.)
    # Then sort, because its much smaller now.
    print(input_bam + " filtering")
    #p1 = subprocess.Popen(["samtools","view","-h", temp_unfiltered_bam ],                                          stdout=subprocess.PIPE)
    #p2 = subprocess.Popen(["grep", "-E", "'(^@)|(XT:Z:)'"],                                      stdin=p1.stdout,  stdout=subprocess.PIPE)
    #p3 = subprocess.Popen(["grep", "-E", "'(^@)|("+corrected_cell_barcode_tag+":Z:)'"],          stdin=p2.stdout,  stdout=subprocess.PIPE)
    #p4 = subprocess.Popen(["samtools", "sort",  "-@", str(threads) ,"-o", peak_anno_bam, "-"],   stdin=p3.stdout,  stdout=subprocess.PIPE)
    #p4.communicate()
    #
    ## ^^^ These commands are not generating the output
    subprocess.check_output("samtools view -h "+temp_unfiltered_bam+" | grep -E '(^@)|(XT:Z:)' | grep -E '(^@)|("+corrected_cell_barcode_tag+":Z:)' | samtools sort  -@ "+str(threads)+" -o "+peak_anno_bam+" -", shell=True)

        
    # check output and cleanup
    if not os.path.exists(peak_anno_bam) or os.stat(peak_anno_bam).st_size == 0 :
        sys.exit("Failed to filter featureCounts output into annotated bam with cmd:\n" +
                  "samtools view -h "+temp_unfiltered_bam+" | grep -E '(^@)|(XT:Z:)' | grep -E '(^@)|("+corrected_cell_barcode_tag+":Z:)' | samtools sort  -@ "+str(threads)+" -o "+peak_anno_bam+" -"   )
    pysam.index(peak_anno_bam)
    if not os.path.exists(peak_anno_bam+".bai")       : sys.exit("No bam index made for file "+peak_anno_bam+" Something went wrong.")       
            
    # Cleanup
    try : 
        os.remove(temp_unfiltered_summary )
        os.remove(temp_unfiltered_bam)
        os.remove(temp_unfiltered)
        os.rmdir(temp_feature_counts_dir)
    except OSError as e:
        warnings.warn("Unable to cleanup after featurecounts. Continuing.")




def count_from_annotated_bams(peak_anno_bams, input_from_dir, counts_root, corrected_cell_barcode_tag, umi_tag) :


    ## Figure out the output.
    counts_files = list()
    if not input_from_dir :
        counts_files.append(counts_root+".tab.gz")
    else :  
        try:  
            os.mkdir(counts_root)
        except OSError:  
            #sys.exit("Couldn't create counts directory (already exists?)"+counts_root)
            print("Counts dir exists already (TEMP)")        
        for peak_anno_bam in peak_anno_bams : 
            sample_base = os.path.splitext(os.path.basename(peak_anno_bam))[0]
            counts_files.append(os.path.join(counts_root, sample_base+".tab.gz"))


    # Process each.
    for n in range(0,len(peak_anno_bams)) : 
        peak_anno_bam = peak_anno_bams[n] 
        counts_file   = counts_files[n]

        print(peak_anno_bam+" counting")
        
        ut_cmd = ["umi_tools", "count",
                "--per-gene", "--gene-tag=XT", "--per-cell", "--extract-umi-method", "tag",
                "--umi-tag="+umi_tag , "--cell-tag="+corrected_cell_barcode_tag,
                "-I", peak_anno_bam,  "-S", counts_file]
        ranok = subprocess.call( ut_cmd )
        if not ranok == 0 :
            sys.exit("Failed to run umi_tools correctly ("+str(ranok)+") with cmd\n"+ " ".join(ut_cmd))

    return(counts_files)






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
    
    if args.polyA_bams and args.peak_anno_bams : # ok, might happen for large dataset.
        sys.exit("Can't specify both --polyA_bams and --peak_anno_bams. If you do have both processed already. " + 
        "Try making the peaks gff file first (if not done alredy)( --polyA_bams with --no_count). "+
        "Then run again with specifying --peaks_gff and --peak_anno_bams")

    if args.peaks_gff and args.polyA_bams :
        sys.exit("Can't specify both --polyA_bams and and --peaks_gff. "+
                 "If there is a peaks gff, polyA bams are no longer needed, use original bams for counting.")

    if args.out_root == "." or args.out_root.endswith("/") :
        sys.exit("--output should be a string, not a directory")



        




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

    try:
        subprocess.call(["umi_tools","--version"])
    except OSError as e:
        sys.exit("Could not find umi_tools in PATH. The umi_tools package should be installed and in PATH.")


    # Version of featureCounts matters, tagging a bam only possible in 1.5.3+
    version_info_str = None
    try: 
        p = subprocess.Popen(["featureCounts","-v"], stderr=subprocess.PIPE)
        version_info_str = p.communicate()[1]
    except OSError as e:
        sys.exit("Could not find featureCounts in PATH. The subread package (v1.5.3 or later) should be installed and in PATH." ) 
    try:    
        version = version_info_str.decode('utf-8').replace('\n','').split()[1].replace('v','')  #Eeeeewwwwww.
        version = version.replace('v','') 
        
        if StrictVersion(version) < StrictVersion('1.5.3') :
                sys.exit("Found featureCounts, but an old version. Need version 1.5.3 or later for this to work.")
    except :
        sys.exit("Coudln't parse the output of 'featureCounts -v'. The subread package (v1.5.3 or later) should be installed and in PATH.")


    print("Tools ok\n")



def quick_bam_check (bam_file, cell_tag, umi_tag, minMAPQ) :
    # is there an index?
    bam_index = bam_file+".bai"
    
    if not os.path.exists(bam_index)       : sys.exit("No bam index "+bam_index+" for file "+bam_file+" Use samtools index input bams.")       

    # Read chr names  (actually, only from bam, so no need.)
    # Look for one of each of those tags in top n reads. Counting won't happen without them!
    top_n = 10000
    n     = 0
    seen_cell_tag = False
    seen_umi_tag  = False
    seen_map_qual_ok = False
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
        sys.exit( "Checked the first "+str(top_n)+" reads of "+bam_file+" and did not see any reads that have min map quality "+str(minMAPQ))


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
