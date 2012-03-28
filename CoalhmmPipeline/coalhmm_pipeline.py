#!/usr/bin/env python

# standard modules
import os, sys, glob, subprocess, socket, math, re, glob, tempfile
from optparse import OptionParser

# temp hack to include Troels code
sys.path.append( "/Users/kasper/projects/ancestral_recomb/coalhmm_pipeline_trunk" )
sys.path.append( "/Users/kasper/projects/ancestral_recomb/" )

# FIXME: hack to override version check for h5
os.putenv("HDF5_DISABLE_VERSION_CHECK", "1")

# ruffus and related tools and utils
from ruffus import *
from RuffusPipeline.RuffusTools import *
from RuffusPipeline.Distribute import PipelineUtils

# # imports for vcf to maf generation
# from CoalhmmPipeline.MAF import MAFIterator
# from CoalhmmPipeline.VCF import SNPIterator
# from CoalhmmPipeline.FastaIndex import FastaIndex
# from CoalhmmPipeline.VCF2MAF import *

# PyTables
from tables import *
from CoalhmmPipeline.OpenHDF import OpenHDF

# imports for input data preparation
from CoalhmmPipeline.Masking import *
from CoalhmmPipeline.Annotation import *
from CoalhmmPipeline.Chunkenizer import Chunkenizer
from CoalhmmPipeline.DummySplitter import DummySplitter
from CoalhmmPipeline.IngroupSplitter import IngroupSplitter
from CoalhmmPipeline.IngroupMergeCriteria import IngroupMergeCriteria
from CoalhmmPipeline.PytablesOutput import PytablesOutput
from CoalhmmPipeline.MafTest import MafTest
from CoalhmmPipeline.MafGroupLengthQualityFilter import MafGroupLengthQualityFilter
from CoalhmmPipeline.MafNContentQualityFilter import MafNContentQualityFilter
from CoalhmmPipeline.CompressCoordinates import compressCoordinates
from CoalhmmPipeline.GenerateLists import *
from CoalhmmPipeline.HDFmerge import mergeHDFfiles

# running coalhmm
from CoalhmmPipeline.CoalHMM import CoalHMM

# processing and storing of coalhmm results
from CoalhmmPipeline.ImportPosteriorTable2 import *
from CoalhmmPipeline.Segmenter import *
from CoalhmmPipeline import Table # FIXME: use something else like what I do in evalutation of simulations in draft_map.py

# utils
from CoalhmmPipeline.PipelineUtils import *

#     root
#      - maps: 
#         main: (alignment) (1 row per segment)
#             colId -key (link to species)
#             alignmentNumber
#             chunk
#             segment
#             score
#             begin (alignment coordinates)
#             end
#         per art:
#             colid -key (link to alignment)
#             chromosome
#             strand
#             begin
#             end
#      - coords (to be compressed)
#         chromosome = UInt16Col(pos=1)
#         strand = Int8Col(pos=2)
#         SpeciesPositionOnPlusStrand = Int64Col(pos=3)
#         alignmentNumber = UInt16Col(pos=4)
#         alignmentPosition = Int64Col(pos=5)

###################################################################################################
### Pipeline control parameters: ##################################################################
###################################################################################################

# list of chromosomes
#CHROMOSOMES = ["chr%d" % x for x in range(1, 23)] + ["chrX"]
CHROMOSOMES = ["chr22"]

# name of reference genome
REFNAME = 'hg18'

# string specifying template alignemnt for each chromosome
TEMPLATEMAFSTRING = '/Users/kasper/NO_BACKUP/bonobo/data/hg18.bonobo_i7.pantro2.ponabe2.mafs/%s.maf'
# ingroup species of template alignment
ORIGINGROUP = ['hg18', 'pantro2', 'bonobo']
# outgroup species of template alignment
ORIGOUTGROUP = ['ponabe2']


REFGENOMEFILE = "DUMMY"
REPLACEMENTS = {}
# REPLACEMENTS = {
#     'hg18':    'human2',
#     'pantro2': 'chimp12',
#     'bonobo':  'bonobo7'
#     }
VCFFILTERS = ['no']

INGROUP = [REPLACEMENTS[k] if k in REPLACEMENTS else k for k in ORIGINGROUP]
OUTGROUP = [REPLACEMENTS[k] if k in REPLACEMENTS else k for k in ORIGOUTGROUP]

MINMAF = 500
MAXGAP = 100
MAXN = 0.7

ANALYSISDIR = "_".join(sorted(INGROUP+OUTGROUP))
DATADIR = "%s/data" % ANALYSISDIR
RESULTDIR = "%s/results" % ANALYSISDIR
LISTDIR = "%s/lists" % ANALYSISDIR
COALHMMDIR = "%s/coalhmm" % ANALYSISDIR
CHUNKDIR = "%s/chunks" % ANALYSISDIR

COALHMMEXEDIR = "./scripts/Simulation/coalhmmOSX"
COALHMMTYPE = "ILS09"
COALHMMOPTIONS = "/Users/kasper/projects/ancestral_recomb/coalhmm.options"
    
MASKANNOTATIONFILE = "/Users/kasper/Dropbox/CoalHMMpipeline/testfiles/anno.txt"
MASK = 'TEST' # Tag defining masking, should encode what file that is used for masking

LISTSIZE = 1e6

###################################################################################################
###################################################################################################
###################################################################################################

def vcf2mafArgumnets():
    """
    Generating function for arguments for maf generation.
    """
    for chrom in CHROMOSOMES:
        yield (TEMPLATEMAFSTRING % chrom, os.path.join(DATADIR, "hg18aln.%s.%s.maf" % (chrom, '_'.join(VCFFILTERS))), REPLACEMENTS, REFNAME, REFGENOMEFILE, VCFFILTERS)

@follows(mkdir(ANALYSISDIR), mkdir(DATADIR))
@files(vcf2mafArgumnets)
def generateMafs(inp, outp, replacements, refName, refGenomeFileName, vcfFilters):
    """
    
    """
    open(outp, 'w').write(open(inp).read())
#     def getVCFfileName(indiv):
#         return indiv + ".vcf"
# 
#     snpSpec = dict()
#     for orig, repl in replacements.items():
#         snpSpec[orig] = SNPIterator(getVCFfileName(repl))
# 
#     vcf2maf = VCF2MAF(MAFIterator(inp), FastaIndex(refGenomeFileName), refName, snpSpec, VCFfilter(vcfFilters))
#     vcf2maf.writeMaf(outp)


@transform(generateMafs, regex(fileRegex), outFile("mask%s" % MASK, outdir=DATADIR), REFNAME, MASKANNOTATIONFILE)
def masking(inp, outp, refName, maskAnnotationFile):
    """
    one chr maf file -> one chr maf file
    """
    refNameRegex = refName + ".*"
    masker = Masking(Annotation(maskAnnotationFile))
    masker.maskMafFile(refNameRegex, inp, outp)


@follows(mkdir(CHUNKDIR))
@transform(masking, regex(fileRegex), outFile("chunkinize", suffix="h5", outdir=CHUNKDIR), CHUNKDIR, INGROUP, OUTGROUP, MINMAF, MAXGAP, MAXN)
def chunkinize(inp, outp, chunkDir, ingroup, outgroup, minMaf, maxGap, maxN):    
    """
    one chr maf file -> one h5 file    
    """
    # This is a bit of a hack: we get the chromosome name from the file name.
    chromName = re.search(r'(chr[\dX]+)', inp).group(1)

    mafQualityFilters = [MafGroupLengthQualityFilter(minMaf, minMaf, ingroup), MafNContentQualityFilter(maxN)]
    chunkQualityFilters = [ChunkNContentQualityFilter(maxN)]

    chunkenizer = Chunkenizer(IngroupSplitter(ingroup, maxGap),
                              IngroupMergeCriteria(maxGap, ingroup),
                              MafTest(ingroup + outgroup),
                              mafQualityFilters,
                              chunkQualityFilters,
                              MafIngroupTruncater(ingroup))

    with OpenHDF(outp, "a") as hdf:
        chunkenizer.chunkenize(inp, PytablesOutput(chunkDir, ingroup, hdf, encodeChrName(chromName)))
    

@follows(mkdir(RESULTDIR))
@collate(chunkinize, regex('.*(h5).*'), CollatedOutputFiles("hdfmerge", "h5", outdir=RESULTDIR), INGROUP, LISTSIZE, CHROMOSOMES)
def buildInputHDF(inp, outp, ingroup, listSize, chromosomes):
    """
    merge alll chunkinize h5 files and compress
    """
    if os.path.exists(outp):
        os.unlink(outp) 
        
    ignore, tmpFileName = tempfile.mkstemp(suffix=".h5")

    # merge
    mergeHDFfiles(inp, tmpFileName)

    # compress coordinates
    hdf = openFile(tmpFileName, 'a')
    for indiv in ingroup:
        compressCoordinates(hdf.root.coordinates[indiv]) 

    # add list info to hdf
    for chrom in chromosomes:
        generateLists(encodeChrName(chrom), listSize, hdf)

    # stop'n'copy garbage collect hdf
    hdf.copyFile(outp) 
    hdf.close()

    os.unlink(tmpFileName)


@follows(mkdir(LISTDIR))
@split(buildInputHDF, regex(fileRegex), outFileGlob("list", suffix="txt", outdir=LISTDIR), LISTDIR, LISTSIZE, CHUNKDIR, CHROMOSOMES)
def writeLists(inp, outp, listDir, listSize, chunkDir, chromosomes):  # list files should be names ...list<chrnr>-<nr>....
    """
    one chunkinize h5 file -> all list files for all chr
    """
    hdf = openFile(inp, 'r')
    outputFileTemplate = getOutputFileName(inp, "list", "_%s_%s", suffix="txt", outdir=LISTDIR)
    for chrom in chromosomes:
        writeAllLists(encodeChrName(chrom), hdf, chunkDir, listDir, outputFileTemplate=outputFileTemplate)
    hdf.close()


@follows(mkdir(COALHMMDIR))
@transform(writeLists, regex(fileRegex), [outFile(COALHMMTYPE, suffix="estimates", outdir=COALHMMDIR), outFile(COALHMMTYPE, suffix="posteriors", outdir=COALHMMDIR)], COALHMMOPTIONS, COALHMMTYPE, COALHMMEXEDIR)
def coalhmm(inp, outp, optionsFileName, coalhmmType, exeDir):
    """
    one list file -> one set of coalhmm result files
    """
#    if distribute(): return

    coalhmm = CoalHMM(optionsFileName)#, coalhmmType, '.')#exeDir)
    commonPrefix = os.path.splitext(outp[0])[0] # get prefix from one of the output file names
    coalhmm.run(inp, commonPrefix)

 
@collate([[coalhmm], buildInputHDF], regex('.*(hg18aln).*'), CollatedOutputFiles("import", "h5", outdir=RESULTDIR, groupDepth=1), REFNAME, CHROMOSOMES)
def importPosteriors(inp, outp, refName, chromosomes):
    """
    all coalhmm result files -> one h5 file with posteriors
    """

    def flatten(lst):
        for elem in lst:
            if type(elem) in (tuple, list):
                for i in flatten(elem):
                    yield i
            else:
                yield elem

    estimateFiles = [i for i in flatten(inp) if i.endswith('.estimates') and os.path.getsize(i)]
    posteriorFiles = [i for i in flatten(inp) if i.endswith('.posteriors') and os.path.getsize(i)]
    inputHDFfileName = [i for i in flatten(inp) if i.endswith('.h5')][0]

    # compute format string for postrior result files:
    commonPrefix = os.path.commonprefix(posteriorFiles)
    commonPrefix =  re.search(r'(.*?_*)[\d_]*$', commonPrefix).group(1)
    commonSuffix = os.path.commonprefix([f[::-1] for f in posteriorFiles])[::-1]
    regex = re.compile('%s(([0-9XY]+)(_+)(\d+))%s' % (commonPrefix, commonSuffix))
    groups = regex.search(posteriorFiles[0]).groups()
    formatStr = commonPrefix + groups[0].replace(groups[1], '%s').replace(groups[3], '%s') + commonSuffix
    
    # create a new hdf not to change mtime of input hdf
    inputHDF = openFile(inputHDFfileName, 'r')
    inputHDF.copyFile(outp) 
    inputHDF.close()

    with OpenHDF(outp, 'a') as hdf:

        for chrom in chromosomes:
            alignmentNumber = encodeChrName(chrom)
            coords = createCoordinatesStructure(hdf.root.coordinates[refName], alignmentNumber)

            for x in hdf.root.lists.where("(alignmentNumber==%s)&(listIndex==1)"%alignmentNumber):
                (listNumber, ignore, ignore, ignore) = x.fetch_all_fields()

                fileName = formatStr % (alignmentNumber, listNumber)
                assert os.path.exists(fileName)
                if os.path.getsize(fileName):
                    importPosterior(hdf, fileName, alignmentNumber, listNumber, coords)

            createSegments(hdf.root.posteriors[chrom])


@files(["test1.in", 'test2.in' ], ["test1.out", 'test2.out'])
def testTask(inp, outp):
    if dpu.distribute(): return
    for i in range(len(inp)):
        open(outp[i], 'w').write(open(inp[i]).read())


@follows(mkdir(RESULTDIR))
@collate(coalhmm, regex('.*(chr).*'), CollatedOutputFiles("tabulate", ["tbl"], outdir=RESULTDIR))
def tabulateEstimates(inp, outp):
    """
    all coalhmm result files for one chr -> one table for one chr
    """
    open(outp, 'w')
    pass
#         userTable = table.Table()
#             #import user.txt anbd params.txt files
#             userfile = "/tmp/results/%i.%i.ILS09_gamma.user.txt"%(alignmentNumber, listNumber)
#             paramsfile = "/tmp/results/%i.%i.ILS09_gamma.params.txt"%(alignmentNumber, listNumber)
#             userTable = userTable.merge(table.Table().load_kv(userfile))
# 
#             userTable.write("/tmp/results/usertable.txt")


if __name__ == "__main__":

    usage = """
    Usage: 
    """

    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-v", "--verbose",
                      dest="verbose",
                      type="int",
                      default=1,
                      help="Level of verbosity")
    parser.add_option("-p", "--printout",
                      action="store_true",
                      default=False,
                      help="Print out what will be run")
    parser.add_option("-g", "--graph",
                      dest="graph",
                      type="string",
                      default=False,
                      help="Task to run")
    parser.add_option("-m", "--multiprocess",
                      dest="multiprocess",
                      type="int",
                      default=1,
                      help="Nr of parallel processes")
    parser.add_option("-t", "--task",
                      dest="task",
                      action="append",
                      type="str",
                      default=[],
                      help="Task to run")
    parser.add_option("-f", "--force",
                      dest="force",
                      action="append",
                      type="str",
                      default=[],
                      help="Task to run")
    parser.add_option("--mode",
                      dest="mode",
                      type="string",
                      default='',
                      help="mode should be 'sge' or 'xgrid'")
    parser.add_option("-l", "--logfile",
                      dest="logfile",
                      type="string",
                      default=None,
                      help="File to log to")
    parser.add_option("-d", "--documentation",
                      dest="documentation",
                      type="string",
                      default=False,
                      help="Print task documentation")
    parser.add_option("-T", "--touch_files_only",
                      action="store_true",
                      default=False,
                      help="Create or update output files only to simulate the running of the pipeline. Does not invoke real task functions to run jobs. This is most useful to force a pipeline to acknowledge that a particular part is now up-to-date. This will not work properly if the identities of some files are not known before hand, and depend on run time. In other words, not recommended if @split or custom parameter generators are being used.")
    (options, args) = parser.parse_args()


    if not options.task and options.force:
        # if no tasks are specified but some forced are we assume that we just want to run
        # the forced ones
        options.task = options.force

    if options.logfile:
        # Set up a specific logger with our desired output level
        my_ruffus_logger = logging.getLogger('My_Ruffus_logger')
        my_ruffus_logger.setLevel(logging.DEBUG)
        # Add the log message handler to the logger
        handler = logging.handlers.RotatingFileHandler(options.logfile)
        logger = my_ruffus_logger.addHandler(handler)
    else:
        logger = stderr_logger

    dpu = PipelineUtils.PipelineUtils(options.mode)
    if options.mode == 'xgrid':
        dpu.setGrid('10.11.101.253', 'birc')
        dpu.addDependencies(glob.glob("../bonobo/coalhmmOSX/*"))
        dpu.addDependencies(["/usr/local/bin/bppseqgen"])
        dpu.addDependencies(glob.glob("pipeline/libbpp*"))

    if options.documentation:
        with open(options.documentation, 'w') as f:
            for n, d in re.findall(r'^@.*?def\s+([^(]*).*?:.*?"""(.*?)"""', open(inspect.stack()[0][1]).read(), re.S | re.M):
                print >>f, n, d
        sys.exit()
        
    if args and options.mode:
        print "Do not try to distribute already distributed tasks..."
        sys.exit()
    elif args:
        # called like this ruffus_pipeline.py task 0,1,1,2 <input_file> <output_file> <extra_option> <extra_option>
        dpu.executeDistributedTask(args)
    else:
        if options.printout:
            pipeline_printout(sys.stderr, options.task, forcedtorun_tasks=options.force, verbose=options.verbose)
        elif options.graph:
            pipeline_printout_graph(options.graph, "jpg", options.task)            
            os.system("qlmanage -p " + options.graph)
        else:
            # called like this ruffus_pipeline.py --mode xgrid
            pipeline_run(options.task, forcedtorun_tasks = options.force,
                         multiprocess=options.multiprocess, verbose=options.verbose,
                         logger=logger, touch_files_only=options.touch_files_only)
