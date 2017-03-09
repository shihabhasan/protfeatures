import time, hashlib

from seq_parameters import *
from pepstat import *
from elemental_sparing import *
from garnier import *
from netcglyc import cMannosylation
from netchop import proteasomalCleavage
from netnglyc import nLinkedGlycosylation
from netphos import *
from prop import arginineLysinePropeptideCleavage
from anchor import bindingRegionsDisordered
from targetp import *
from tmhmm import transmembraneHelices
from signalp import signalPeptides
from bepipred import linearBcellEpitopes
from immuno import classIimmunogenicity


def features(seqID, sequence, features_list):
    parameters=""   
    seq=sequence.upper().replace("X","").replace("*","").replace("U","C").replace("B","D")
    inFileName=seqID.replace("|","_").replace(":","_")+"_"+hashlib.md5(time.asctime()).hexdigest()
    inFasta=open(inFileName, "w")
    inFasta.write(">"+seqID.replace("|","_")+"\n"+seq+"\n")
    inFasta.close()
    for feature in features_list:
        feature=eval(feature)
        parameters=parameters+feature(inFileName)+","
    
    os.remove(inFileName)
    parameters = parameters.rstrip(",")
    
    return parameters

    
