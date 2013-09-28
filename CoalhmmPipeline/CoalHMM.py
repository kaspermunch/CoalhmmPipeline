
import os, sys
import subprocess

# method = 'ILS09_gamma'
# listFileName = "/tmp/lists/tmp.lst" # name of list file
# resultFilePrefix = "/tmp/results/22.1" # base path of all result files (containing list id)
# coalHMMdir="/users/tth/local/coalhmm/CoalHMM/"

class CoalHMM(object):

    def __init__(self, executable, optionsFileName, ingroup, outgroup, **parameters):#, coalHMMdir):
        self.executable = executable
        self.optionsFileName = optionsFileName
        self.parameters = parameters
#        self.coalHMMdir = coalHMMdir
        self.success = True

        self.speciesSpec = 'species1=%s species2=%s species3=%s outgroup=%s' % (ingroup[0], ingroup[1], ingroup[2], outgroup[0])

    def run(self, listFileName, outputFilePrefix):

        extraArgs = ' '.join(['='.join(map(str, t)) for t in self.parameters.items()])

#         # add executable dir name to LD path in case library files are there:
#         os.environ['LD_LIBRARY_PATH'] += ':' + os.path.dirname(self.executable)

        cmd = "%s --noninteractive=yes param=OPTIONSFILE %s input.sequence.list=LISTFILE input.sequence.multiparts.prefix= input.sequence.multiparts.reset=yes optimization.profiler=BASENAME.profiler optimization.message_handler=BASENAME.messages output.posterior.states=none output.hidden_states=BASENAME.states output.hidden_states.divergences=BASENAME.divergences output.posterior.values=BASENAME.posteriors output.estimated.parameters=BASENAME.params output.userfriendly.parameters=BASENAME.estimates %s" % (self.executable, self.speciesSpec, extraArgs)

        cmd = cmd.replace('BASENAME', outputFilePrefix).replace('LISTFILE', os.path.abspath(listFileName)).replace('OPTIONSFILE', os.path.abspath(self.optionsFileName))

        #p = subprocess.Popen(cmd, env=os.environ, shell=True, cwd=self.coalHMMdir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        p = subprocess.Popen(cmd, env=os.environ, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()

        if p.returncode != 0:
            self.success = False

        errorList = []

        if any(error in stdout+stderr for error in errorList):
            self.success = False

        with open("%s.out" % outputFilePrefix, 'w') as f:
            f.write(stdout) 
        with open("%s.err" % outputFilePrefix, 'w') as f:
            f.write(stderr) 

    def cleanUp(self):

        # parse files




        # cleanup
        if not self.success:
            print "Faliure"
            for o in glob.glob(self.optionsFileName + '*'):
                os.remove(o)

        







