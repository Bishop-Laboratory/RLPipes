#!/usr/bin/env python
# Script Name: QmRLFS-finder.py (version v2)
# Description: QmRLFS-finder.py is program to predict 
#              R-loop forming regions in any DNA sequences 
#              (in FASTA file)
# Created By:  Piroon Jenjaroenpun and Thidathip Wongsurawat
# e-mail:      piroonj@gmail.com
# Update:      29-Aug-2019
################################
##History
## QmRLFS-finder v2.0
##   Fix: program for Python2 and Python3
## QmRLFS-finder v1.5
##   Fix: coordinate name in fasta result and name in bed result
## QmRLFS-finder v1.4
##   add ENSEMBL format coordinate recognition
##   add UCSC format coordinate recognition
##   add Bedtools format coordinate recognition
##------------------------------
## QmRLFS-finder v1.3
##   increase speed of finding REZ
##   introduce pre-filtering non RLFS sequence
##------------------------------
## QmRLFS-finder v1.2
##  add --log option
##  add --strand option
################################

import datetime, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random, re, string, sys, argparse, math, os
import textwrap

class GenerateWindows():
    def __init__(self, length=1000, step=10):
        self.length = length
        self.step = step
    def __iter__(self): # must return an iterator!
        return iter(range(0,self.length,self.step))

prog_name = os.path.basename(sys.argv[0])
version = "v2.0"
today = datetime.datetime.today()
parser = argparse.ArgumentParser(prog="python {}".format(prog_name), \
                                formatter_class=argparse.RawDescriptionHelpFormatter,\
                                description=textwrap.dedent('''\
                                    Description: 
                                     {} (version {}) is program to predict R-loop 
                                     forming regions in any DNA sequences'''.format(prog_name, version)), \
                                usage=textwrap.dedent('''%(prog)s [options] -i <fasta> -o <fname>''')
                                )

in_params = parser.add_argument_group(' ------- Input file -----------')
in_params.add_argument('-i','--infile', dest='infile', metavar='<fasta>',\
                       type=str, default="", \
                       help="input in FASTA file (default: stdin)")

model_params = parser.add_argument_group(' ------- QmRLFS Model ---------', 'To set QmRLFS model(s) for RLFS finding')
model_params.add_argument('-m', '--model', type=str, action='store',
                    dest='models', default="m1,m2", metavar='<str>',
                    help='choose QmRLFS model(s) for RLFS finding. The user chooses single model m1 or m2 or multiple models e.g. m1,m2. The detailed information is available on-line at http://rloop.bii.a-star.edu.sg/?pg=program. (default: %(default)s)') 

output_params = parser.add_argument_group(' ------- Output format --------','To define output name and formats')
output_params.add_argument('-o','--outfile', dest='outfile', metavar='<fname>',\
                            type=str, default="", \
                            help="write RLFSs to file <fname> (default: stdout)")
output_params.add_argument('-table', action='store_true', dest='table', \
						   help='create RLFS table file format (selected by default)')
output_params.add_argument('-bed', action='store_true', dest='bed', \
						   help='create BED file format')
output_params.add_argument('-fa', action='store_true', dest='fa', \
						   help='create FASTA file format')
output_params.add_argument('-custom_track', action='store_true', dest='custom_track', \
						   help='create CUSTOM TRACK file format for visualizing on UCSC browser')
                           
other_params = parser.add_argument_group(' ------- Others ---------------','The options for adjusting output')
other_params.add_argument('--strand', type=str, dest='strand', choices=['+', '-'], default="+", \
						   help='genomic strand (default: +)')
other_params.add_argument('--coordinate', type=str, dest='coordinate', metavar='<chr:start-end>', \
						   help='re-calculate RLFS start and end using provided genomic position e.g. chr1:10000-11000 (chr:start-end)')
other_params.add_argument('--coordinate_name', action='store_true', dest='coordinate_name', \
						   help='re-calculate RLFS start and end using genomic position in FASTA name. the FASTA name should contains "chr:start-end" e.g. chr1:10000-11000 or "chromosome:assembly:chr:start:end" e.g. chromosome:GRCh38:1:10000:11000')
other_params.add_argument('--coordinate_base_zero', action='store_true', dest='coordinate_base_zero', \
						   help='The first base in a chromosome is numbered 0.')
other_params.add_argument('--track_name', type=str, dest='track_name',metavar='<str>', default='RLFS', \
						   help='define track name  for visualizing on UCSC browser. (default: %(default)s)')
other_params.add_argument('--track_description', type=str, dest='track_description',metavar='<str>', \
                           default='RLFS prediction via QmRLFS-finder', \
						   help='define track description for visualizing on UCSC browser. (default: %(default)s)')
other_params.add_argument('--keep_name', action='store_true', dest='keep_name', \
						   help='keep sequence name in output file, when use --coordinate or --coordinate_name.')
other_params.add_argument('--quick', action='store_true', dest='speed', \
                           help='increase speed of searching REZ by skipping REZ length maximization. In general, RLFS regions searched by this option are shorter than default.')
other_params.add_argument('--log', action='store_true', dest='log', \
						   help='Write stdout in logfile')
other_params.add_argument('--verbose', action='store_true', dest='verbose', \
						   help='Print verbose output')
other_params.add_argument('--version', '-v', action='store_true', dest='version', \
						   help='Print version')

args = parser.parse_args()

## Manage input and parameters ############
if args.version:
	exit("QmRLFS-finder version {}".format(version))

if args.infile == "":
    infile = sys.stdin
else:
    infile = open(args.infile,"rU")
if infile.name == "<stdin>":
    try:
        infile.seek(0)
    except:
        print("""\n-----------------------------------\n !!! There is no input File !!!\n-----------------------------------\n""")
        exit(parser.print_help())

table_out = False
if args.outfile == "":
    if len(re.findall("True", repr([str(args.bed), str(args.fa), str(args.custom_track), str(args.table)]))) > 1:
        exit("\nIf -o <fname> or --outfile <fname> is not defined, please select only \none output format of -table, -fa, -bed, or -custom_track")
    if not args.table and not args.fa and not args.bed and not args.custom_track:
        table_out = True
    if args.table or table_out:
        outfile_table = sys.stdout
    elif args.fa:
        outfile_fa = sys.stdout
    elif args.bed:
        outfile_bed = sys.stdout
    elif args.custom_track:
        outfile_custom_track = sys.stdout
else:
    if not args.table and not args.fa and not args.bed and not args.custom_track:
        table_out = True
    if args.table or table_out:
        outfile_table = open("{}.out.table.txt".format(args.outfile),"w")
    if args.fa:
        outfile_fa = open("{}.out.fa".format(args.outfile),"w")
    if args.bed:
        outfile_bed = open("{}.out.bed".format(args.outfile),"w")
    if args.custom_track:
        outfile_custom_track = open("{}.out.customTrack.bed".format(args.outfile),"w")
supported_models = set(["m1","m2"])
input_models = set(args.models.split(","))
unknown_models = set(input_models.difference(supported_models))
if len(unknown_models)>0:
    exit('\n argument -m or --model does not support "{}"'.format(",".join(list(unknown_models))))
    
########################################
## general setting
dateformat = "%a %b %d %Y %H:%M:%S "
re_position = re.compile(r"(\d+)-(\d+)")
re_pyr = re.compile(r"^G")
re_coordinate_ucsc = re.compile(r"(\w+):(\d+)-(\d+).*strand.([+-])")
re_coordinate_bedtools = re.compile(r"(\w+):(\d+)-(\d+)[\.\(]([+-])[\.\)]")
re_coordinate_no_strand = re.compile(r"(\w+):(\d+)-(\d+)")
re_coordinate_ensembl = re.compile(r"chromosome:\w+:(\w+):(\d+):(\d+):(-*\d+)")
re_coordinate_ensembl_no_strand = re.compile(r"chromosome:\w+:(\w+):(\d+):(\d+)")
########################################
## parameters
models = {}
if "m1" in input_models:
    models['m1'] = re.compile(r"G{3,}[ATCGU]{1,10}?G{3,}(?:[ATCGU]{1,10}?G{3,}){1,}?")
if "m2" in input_models:
    models['m2'] = re.compile(r"G{4,}(?:[ATCGU]{1,10}?G{4,}){1,}?")
params = {'regex_models':models,
#---------RIZ parameters---------
          'min_perc_g_riz':50, 
#---------Linker parameters---------
          'num_linker':50, 
#---------REZ parameters---------
          'window_step': 100, 
          'max_length_rez':2000, 
          'min_perc_g_rez':40
          }
########################################
# HEAD line of RLFS format
tableheaders_combined = ["#Name", "model", "location", "start_RIZ", "end_RIZ", "length_RIZ", "G_RIZ",
						 "3Gs_RIZ", "4Gs_RIZ", "perc_G_RIZ", "sequence_RIZ",
						 "Linker", "start_REZ", "end_REZ", "length_REZ", "G_REZ",
						 "3Gs_REZ", "4Gs_REZ", "perc_G_REZ", "sequence_REZ", "strand"]
########################################################
def percent_g(seq): # checked
    """Returns the percentage of G bases in sequence seq"""
    return float("%.2f"%((seq.count("G") / float(len(seq))) * 100))
########################################################
def riz_search(seq, model):
    """Searches and returns RLFS RIZ regions in the given sequence seq"""
    result_list = []
    for result in params['regex_models'][model].finditer(seq):
        start = result.start()
        end = result.end()
        riz_seq = result.group(0)
        perc_g = percent_g(riz_seq)
        if perc_g >= params['min_perc_g_riz']:
            dict_result = {"start": start,
                           "end": end,
                           "length": end - start,
                           "G": riz_seq.count("G"),
                           "G3s": riz_seq.count("GGG"),
                           "G4s": riz_seq.count("GGGG"),
                           "perc_g": perc_g,
                           "model" : model,
                           "seq": riz_seq}
            result_list.append(dict_result)
    return result_list
########################################################
def rez_search(riz_end, seq, seqlength): 
    """Searches and returns RLFS RIZ regions in the given sequence seq"""
    rez_regions = {}
    start_seq = 0
    end_seq = 0
    max_length = 0
    preFilteringRLFS = 0
    int_step = 10
    int_window = 100
    ## to pre-filter RLFS before search ------------------------
    if args.speed:
        window_start = GenerateWindows(len(seq[:1000]), int_step)
        # print len(seq[:2100]), seq
        for i_win in window_start:
            tmpSeq = seq[i_win:i_win+int_window]
            perc_tmpSeq = percent_g(tmpSeq)
            # print i_win, i_win+int_window, " : ", perc_tmpSeq
            if perc_tmpSeq >= params['min_perc_g_rez']:
                preFilteringRLFS = 1
        # print "-----------------------------------"
        if preFilteringRLFS == 0:
            window_start = GenerateWindows(len(seq[:150]), 1)
            for i_win in window_start:
                tmpSeq = seq[i_win:i_win+int_window]
                perc_tmpSeq = percent_g(tmpSeq)
                # print i_win, i_win+int_window, " : ", perc_tmpSeq
                if perc_tmpSeq >= 35:
                    preFilteringRLFS = 1
            # print "-----------------------------------"
            if preFilteringRLFS == 0:
                # print "not pass"
                return {}
    ##------------------------  
    # print "pass pre-filter"
    for i in range(params['num_linker']):
        seq_frag = str(seq[i:i + params['window_step']])
        if seq_frag:
            if args.speed:
                rez_loop_search = params['max_length_rez']/2
            else:
                rez_loop_search = params['max_length_rez']
            for j in range(i + params['window_step'], i + rez_loop_search + 1):
                seq_tmp = seq[i:j]
                if seq_tmp and percent_g(seq_tmp) >= params['min_perc_g_rez'] and seqlength > j \
                            + riz_end and max_length < len(seq_tmp):
                        max_length = len(seq_tmp)
                        start_seq = i
                        end_seq = j
        if args.speed and max_length != 0:
            break
    if max_length:
        rez_regions = get_rez_element(riz_end, start_seq, end_seq, seq)
    return rez_regions
########################################################
def get_rez_element(riz_end, rez_start, rez_end, seq):
    """Returns a dict with the elements of the specified REZ region."""
    rez_seq = seq[rez_start:rez_end]
    return {"start": riz_end + rez_start,
            "end": riz_end + rez_end,
            "length": rez_end - rez_start,
            "G": rez_seq.count("G"),
            "G3s": rez_seq.count("GGG"),
            "G4s": rez_seq.count("GGGG"),
            "perc_g": percent_g(rez_seq),
            "seq": rez_seq}
########################################################
def output_rlfs(riz, rez, seq, strand, desc, model):
    """Takes the RIZ and REZ regions from an RLFS as input and generates
    the necessary output.
    """
    global outfile_table, outfile_fa, outfile_bed, outfile_custom_track, table_out
    seqlength = len(seq)
    riz_rez_sorted = sorted([riz['start'], riz['end'], rez['start'], rez['end']])
    linker_length = riz_rez_sorted[2] - riz_rez_sorted[1]
    linker_seq = seq[riz_rez_sorted[1]:riz_rez_sorted[2]]
        # no need to adjust "+"
        #  RIZ start   RIZ end  REZ start                     REZ end
        #       |         |         |                            |
        # 5' ------------------------------------------------------ 3'
    if strand == "-":
        # to adjust the start and end coordinates applied to the
        # reverse complement sequence.  
        #  REZ start                   REZ end  RIZ start    RIZ end
        #       |                         |         |           |
        # 3' ------------------------------------------------------ 5'
        RIZ_start = seqlength - riz["end"]
        RIZ_end = seqlength - riz["start"]
        REZ_start = seqlength - rez["end"]
        REZ_end = seqlength - rez["start"]
        riz["start"] = RIZ_start
        riz["end"] = RIZ_end
        rez["start"] = REZ_start
        rez["end"] = REZ_end
    
    ## get OUTPUT variables
    ## get chrom_i, start_i, and name
    strand_i = strand
    if args.coordinate_name or args.coordinate:
        ## Extract chromosome and start
        if args.coordinate_name:
            m_ucsc_coord_ucsc = re_coordinate_ucsc.search(desc.replace(",",""))
            m_ucsc_coord_bedtools = re_coordinate_bedtools.search(desc.replace(",",""))
            m_ucsc_coord_no_strand = re_coordinate_no_strand.search(desc.replace(",",""))
            m_ensembl_coord = re_coordinate_ensembl.search(desc.replace(",",""))
            m_ensembl_coord_no_strand = re_coordinate_ensembl_no_strand.search(desc.replace(",",""))
            if m_ucsc_coord_ucsc:
                (chrom_i, start_i, end_i, strand_i) = m_ucsc_coord_ucsc.groups()
            elif m_ucsc_coord_bedtools:
                (chrom_i, start_i, end_i, strand_i) = m_ucsc_coord_bedtools.groups()
            elif m_ucsc_coord_no_strand:
                (chrom_i, start_i, end_i) = m_ucsc_coord_no_strand.groups()
            elif m_ensembl_coord:
                (chrom_i, start_i, end_i, strand_i) = m_ensembl_coord.groups()
            elif m_ensembl_coord_no_strand:
                (chrom_i, start_i, end_i) = m_ensembl_coord_no_strand.groups()
            else:
                print( "\n--coordinate_name, the FASTA name should contains \n \"chr:start-end\" e.g. chr1:10000-11000 \n \"chromosome:assembly:chr:start:end\" e.g. chromosome:GRCh38:1:10000:11000\n UCSC format -> \"chr:start-end strand=+\" e.g. chr1:10000-11000 strand=+\n Bedtools format -> \"chr:start-end(strand)\" e.g. chr1:10000-11000(+)\n ENSEMBL format -> \"chromosome:assembly:chr:start:end:strand\" e.g. chromosome:GRCh38:1:10000:11000:1")
                exit()
        else:
            try:
                (chrom_i, start_i, end_i) = re_coordinate_no_strand.search(args.coordinate.replace(",","")).groups()
            except: 
                exit("\n--coordinate should be in chr:start-end e.g. chr1:10000-11000")
        start_i = int(start_i) 
        end_i = int(end_i)
        if args.coordinate_base_zero:
            start_i = start_i+1
        if end_i < start_i:
            exit("\nIn --coordinate chr:start-end and --coordinate_name, start must less than end.")
        elif seqlength != end_i-start_i+1 and args.coordinate_base_zero:
            exit("In --coordinate chr:start-end and --coordinate_name, end-start+1 ({}-{}+1={}) must equal to input sequence length ({})\nexclude --coordinate_base_zero option may solve the problem.".format(end_i,start_i,end_i-start_i+1,seqlength))
        elif seqlength != end_i-start_i+1:
            exit("\nIn --coordinate chr:start-end and --coordinate_name, end-start+1 ({}-{}+1={}) must equal to input sequence length ({})\ninclude --coordinate_base_zero option may solve the problem.".format(end_i,start_i,end_i-start_i+1,seqlength))
        start_i = start_i-1
    else:
        chrom_i = desc
        start_i = 0
    ## Set the output strand specific values
    if strand == "+":
        start = start_i + int(riz['start'])
        end = start_i + int(rez['end'])
        blockSize = "{},{}".format(riz['length'],rez['length'])
        blockStarts = "{},{}".format(0,int(rez['start'])-int(riz['start']))
    else:
        start = start_i + int(rez['start'])
        end = start_i + int(riz['end'])
        blockSize = "{},{}".format(rez['length'],riz['length'])
        blockStarts = "{},{}".format(0,int(riz['start'])-int(rez['start']))
        
    if args.keep_name:
        name = desc
    else:
        name = "{}.{}.{}".format(chrom_i,end-start,start)
        
    if args.coordinate_name or args.coordinate:
        location = "{}:{}-{}".format(chrom_i,start+1,end)
        fasta_comment = "{}.{} ({}:{}-{}) (model={};strand={})".format(name,riz['model'],chrom_i,start+1,end, model, strand)
    else:
        location = "{}:{}-{}".format(desc,start+1,end)
        fasta_comment = "{}.{} (model={};strand={})".format(name,riz['model'], model, strand)
    ## Set the output fix value
    score = itemRgb = "0"
    strand = strand
    blockCount = "2"
    ## Write output
    if args.fa:
        fa_records.append(SeqRecord(Seq(seq[min(riz_rez_sorted):max(riz_rez_sorted)]), id=fasta_comment, description=""))
    if args.bed:
        name = name.replace(" ",".")
        chrom_i = chrom_i.replace(" ",".")
        outfile_bed.write("{}\n".format("\t".join(map(str, [chrom_i, start, end, "{}.{}".format(name,riz['model']), score, strand, start, end, itemRgb, blockCount, blockSize, blockStarts]))))
    if args.custom_track:
        name = name.replace(" ",".")
        chrom_i = chrom_i.replace(" ",".")
        outfile_custom_track.write("{}\n".format("\t".join(map(str, [chrom_i, start, end, "{}.{}".format(name,riz['model']), score, strand, start, end, itemRgb, blockCount, blockSize, blockStarts]))))

    if args.table or table_out:
            outfile_table.write("{}\n".format("\t".join(map(str, [desc, riz['model'], location, riz["start"], riz["end"], riz["length"],
                    riz["G"], riz["G3s"], riz["G4s"], riz["perc_g"],
                    riz["seq"], linker_length, rez["start"],
                    rez["end"], rez["length"], rez["G"], rez["G3s"],
                    rez["G4s"], rez["perc_g"], rez["seq"], strand]))))
########################################################
def search_rlfs(model, seq, strand, desc): 
    """Searches the specified sequence seq for RLFS, and prints the information
       to a file in CSV format.
    """
    global outfile_table, outfile_fa, outfile_bed, outfile_custom_track
    global head, time_p, outlog
    time_process = "#==0==0==0==0==0==0==0==0==0==#"
    ## find RIZ regions
    riz_list = riz_search(seq, model)
    seqlength = len(seq)
    number_rlfs = 0
    if riz_list:  
        for riz in riz_list:
	    ## find REZ for each riz
            if args.outfile != "" and args.verbose:
                time_point = int(float(riz['start'])/seqlength*100/3.3)
                if not time_point in set(time_p) and time_point>0:
                    if args.log:
                        outlog.write( time_process[time_p[-1]:time_point] )
                        outlog.flush()
                    else:
                        sys.stdout.write( time_process[time_p[-1]:time_point] )
                        sys.stdout.flush()
                    time_p.append(time_point)
            rez = rez_search(riz["end"], seq[riz["end"]:], seqlength)
            if rez:
                ## print head file
                if head:
                    if args.table or table_out:
                        outfile_table.write("{}\n".format("\t".join(tableheaders_combined)))
                    if args.custom_track:
                        outfile_custom_track.write('track name="{}" description="{}" visibility=4 colorByStrand="5,190,255 244,0,224"\n'.format(args.track_name, args.track_description))
                    head = False
                
                ## RLFS report
                output_rlfs(riz, rez, seq, strand, desc, model)
                number_rlfs += 1

        if args.outfile != "" and args.verbose:
            if args.log:
                outlog.write( time_process[time_p[-1]:] )
            else:
                print( time_process[time_p[-1]:])
    else:
        if args.outfile != "" and args.verbose:
            if args.log:
                outlog.write( time_process[time_p[-1]:] )
            else:
                print( time_process[time_p[-1]:])
     
    if args.outfile != "" and args.verbose:
        if args.log:
            outlog.write( "\n{}> Found {} RLFSs\n".format(" "*6,number_rlfs) )
        else:
            print( "\n{}> Found {} RLFSs".format(" "*6,number_rlfs))

########################################################
#                  MAIN
########################################################
def main():
    global outfile_table, outfile_fa, outfile_bed, outfile_custom_track
    global fa_records, head, time_p, outlog
    # iterate over all fasta sequences in the input file
    if args.outfile != "":
        if args.log:
            outlog = open("{}.log".format(args.outfile), "wb")
        if args.log:
            outlog.write( "{} (version {})\n {} {}\n".format( prog_name, version, "run on", today.strftime(dateformat) ) )
            outlog.write( " command line: python {}\n".format(" ".join(sys.argv) ))
        else:
            print( "{} (version {})\n {} {}".format( prog_name, version, "run on", today.strftime(dateformat) ))
            print (" command line: python {}".format(" ".join(sys.argv)))

    head = True
    fa_records = []
    for record in SeqIO.parse(infile, "fasta"):
        if args.outfile != "" and args.verbose:
            if args.log:
                outlog.write( "\n{}Search RLFS in {}\n".format(">"*2,record.description) )
            else:
                print ("\n{}Search RLFS in {}".format(">"*2,record.description))
                
        if args.strand:
            genomic_strand = args.strand
        else:
            genomic_strand = "+"
            
        if args.coordinate_name:
            m_ucsc_coord_ucsc = re_coordinate_ucsc.search(record.description.replace(",",""))
            m_ucsc_coord_bedtools = re_coordinate_bedtools.search(record.description.replace(",",""))
            m_ensembl_coord = re_coordinate_ensembl.search(record.description.replace(",",""))
            if m_ucsc_coord_ucsc:
                genomic_strand = m_ucsc_coord_ucsc.group(4)
            elif m_ucsc_coord_bedtools:
                genomic_strand = m_ucsc_coord_bedtools.group(4)
            elif m_ensembl_coord:
                if m_ensembl_coord.group(4) == "1":
                    genomic_strand = "+"
                elif m_ensembl_coord.group(4) == "-1":
                    genomic_strand = "-"
                else:
                    print ('\nEnsemble strand format should be 1 or -1')
                    exit()
                
        for model in params["regex_models"]:
            #### positive strand
            time_p = [1]
            tmp_start_time = time.time()
            if args.outfile != "" and args.verbose:
                if args.log:
                    outlog.write( "\n{}Model: {}\n\n".format(" "*4,model) )
                    outlog.write( "{}- Search RLFS in positive strand\n\n{}> Process\n".format(" "*4," "*6) )
                    outlog.write( "{}#--1--2--3--4--5--6--7--8--9--# %\n".format(" "*8) )
                    outlog.write( "{}#".format(" "*8))
                    outlog.flush()
                else:
                    print ("\n{}Model: {}\n".format(" "*4,model))
                    print ("{}- Search RLFS in positive strand\n\n{}> Process".format(" "*4," "*6))
                    print ("{}#--1--2--3--4--5--6--7--8--9--# %".format(" "*8))
                    sys.stdout.write( "{}#".format(" "*8))
            #
            if genomic_strand == '-':
                seq = str(record.seq.reverse_complement().upper())
            else:
                seq = str(record.seq.upper())
            search_rlfs(model, seq, "+", record.description)
            #
            if args.outfile != "" and args.verbose:
                if args.log:
                    outlog.write(  "\n{0}> Time used: {1:.2f} mins\n\n".format(" "*6,(time.time() - tmp_start_time) / 60) )
                else:
                    print ("\n{0}> Time used: {1:.2f} mins\n".format(" "*6,(time.time() - tmp_start_time) / 60))
            
            #### negative strand
            time_p = [1]
            tmp_start_time = time.time()
            if args.outfile != "" and args.verbose:
                if args.log:
                    outlog.write( "{}- Search RLFS in negative strand\n\n{}> Process\n".format(" "*4," "*6) )
                    outlog.write( "{}#--1--2--3--4--5--6--7--8--9--# %\n".format(" "*8) )
                    outlog.write( "{}#".format(" "*8) ) 
                    outlog.flush()
                else:
                    print ("{}- Search RLFS in negative strand\n\n{}> Process".format(" "*4," "*6))
                    print ("{}#--1--2--3--4--5--6--7--8--9--# %".format(" "*8))
                    sys.stdout.write( "{}#".format(" "*8) )
            #
            if genomic_strand == '-':
                seq = str(record.seq.upper())
            else:
                seq = str(record.seq.reverse_complement().upper())
            search_rlfs(model, seq, "-", record.description)
            #
            if args.outfile != "" and args.verbose:
                if args.log:
                    outlog.write(  "\n{0}> Time used: {1:.2f} mins\n\n".format(" "*6,(time.time() - tmp_start_time) / 60) )
                else:
                    print ("\n{0}> Time used: {1:.2f} mins\n".format(" "*6,(time.time() - tmp_start_time) / 60))
    if args.fa:
        SeqIO.write(fa_records, outfile_fa, "fasta")

if __name__ == "__main__":
    start_time = time.time()
    main()
    infile.close()
    if args.table or table_out:
        outfile_table.close()
    if args.bed:
        outfile_bed.close()
    if args.custom_track:
        outfile_custom_track.close()
    if args.log:
        outlog.close()
    if args.outfile != "":
        print ("\nTime used: {0:.2f} mins".format((time.time() - start_time) / 60))
