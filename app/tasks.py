from __future__ import absolute_import
from celery import task

import os, sys, subprocess
import time
import sqlite3
import hashlib
sys.path.append('/home/shihab/protfeatures/scripts')
from features import features
from duplicate_seq_remover import duplicate_seq_remover
from Bio import SeqIO
from StringIO import StringIO
from sklearn import svm, preprocessing
import numpy as np
from datetime import datetime


featuresDic={'seqLength':'Length of sequence',
'molecularWeight':'Molecular Weight',
'averageResidueWeight':'Average Residue Weight',
'alanineCount':'Count of alanine',
'cysteineCount':'Count of cysteine',
'asparticAcidCount':'Count of aspartic acid',
'glutamicAcidCount':'Count of glutamic acid',
'phenylalanineCount':'Count of phenylalanine',
'glycineCount':'Count of glycine',
'histidineCount':'Count of histidine',
'isoleucineCount':'Count of isoleucine',
'lysineCount':'Count of lysine',
'leucineCount':'Count of leucine',
'methionineCount':'Count of methionine',
'asparagineCount':'Count of asparagine',
'prolineCount':'Count of proline',
'glutamineCount':'Count of glutamine',
'arginineCount':'Count of arginine',
'serineCount':'Count of serine',
'threonineCount':'Count of threonine',
'valineCount':'Count of valine',
'tryptophanCount':'Count of tryptophan',
'tyrosineCount':'Count of tyrosine',
'tinyMoleCount':'Count of tiny mole',
'smallMoleCount':'Count of small mole',
'alaninePercentage':'Percentage of alanine',
'cysteinePercentage':'Percentage of cysteine',
'asparticAcidPercentage':'Percentage of aspartic acid',
'glutamicAcidPercentage':'Percentage of glutamic acid',
'phenylalaninePercentage':'Percentage of phenylalanine',
'glycinePercentage':'Percentage of glycine',
'histidinePercentage':'Percentage of histidine',
'isoleucinePercentage':'Percentage of isoleucine',
'lysinePercentage':'Percentage of lysine',
'leucinePercentage':'Percentage of leucine',
'methioninePercentage':'Percentage of methionine',
'asparaginePercentage':'Percentage of asparagine',
'prolinePercentage':'Percentage of proline',
'glutaminePercentage':'Percentage of glutamine',
'argininePercentage':'Percentage of arginine',
'serinePercentage':'Percentage of serine',
'threoninePercentage':'Percentage of threonine',
'valinePercentage':'Percentage of valine',
'tryptophanPercentage':'Percentage of tryptophan',
'tyrosinePercentage':'Percentage of tyrosine',
'tinyMolePercentage':'Percentage of tiny mole',
'smallMolePercentage':'Percentage of small mole',
'alanineDayhoffStat':'Dayhoff statistic of alanine',
'cysteineDayhoffStat':'Dayhoff statistic of cysteine',
'asparticAcidDayhoffStat':'Dayhoff statistic of aspartic acid',
'glutamicAcidDayhoffStat':'Dayhoff statistic of glutamic acid',
'phenylalanineDayhoffStat':'Dayhoff statistic of phenylalanine',
'glycineDayhoffStat':'Dayhoff statistic of glycine',
'histidineDayhoffStat':'Dayhoff statistic of histidine',
'isoleucineDayhoffStat':'Dayhoff statistic of isoleucine',
'lysineDayhoffStat':'Dayhoff statistic of lysine',
'leucineDayhoffStat':'Dayhoff statistic of leucine',
'methionineDayhoffStat':'Dayhoff statistic of methionine',
'asparagineDayhoffStat':'Dayhoff statistic of asparagine',
'prolineDayhoffStat':'Dayhoff statistic of proline',
'glutamineDayhoffStat':'Dayhoff statistic of glutamine',
'arginineDayhoffStat':'Dayhoff statistic of arginine',
'serineDayhoffStat':'Dayhoff statistic of serine',
'threonineDayhoffStat':'Dayhoff statistic of threonine',
'valineDayhoffStat':'Dayhoff statistic of valine',
'tryptophanDayhoffStat':'Dayhoff statistic of tryptophan',
'tyrosineDayhoffStat':'Dayhoff statistic of tyrosine',
'carbonCount':'Count of carbon sparing',
'nitrogenCount':'Count of nitrogen sparing',
'sulphurCount':'Count of sulphur sparing',
'oxygenCount':'Count of oxygen sparing',
'hydrogenCount':'Count of hydrogen sparing',
'carbonSparingAverage':'Average carbon sparing',
'nitrogenSparingAverage':'Average nitrogen sparing',
'sulphurSparingAverage':'Average sulphur sparing',
'oxygenSparingAverage':'Average oxygen sparing',
'hydrogenSparingAverage':'Average hydrogen sparing',
'aromaticity':'Aromaticity',
'instabilityIndex':'Instability Index',
'isoelectricPoint':'Isoelectric Point',
'gravy':'Grand average of hydropathy (GRAVY)',
'charge':'Charge',
'molarExtinctionCoefficient':'Molar Extinction Coefficient A280',
'absobance':'Absobance A280',
'probabilityOfExpression':'Probability of Expression Inclusion Bodies',
'aliphaticMoleCount':'Count of aliphatic mole',
'aromaticMoleCount':'Count of aromatic mole',
'polarMoleCount':'Count of polar mole',
'nonPolarMoleCount':'Count of non polar mole',
'chargedMoleCount':'Count of charged mole',
'acidicMoleCount':'Count of acidic mole',
'basicMoleCount':'Count of basic mole',
'aliphaticMolePercentage':'Percentage of aliphatic mole',
'aromaticMolePercentage':'Percentage of aromatic mole',
'polarMolePercentage':'Percentage of polar mole',
'nonPolarMolePercentage':'Percentage of non polar mole',
'chargedMolePercentage':'Percentage of charged mole',
'acidicMolePercentage':'Percentage of acidic mole',
'basicMolePercentage':'Percentage of basic mole',
'helixFraction':'Secondary helix fraction',
'turnFraction':'Secondary turn fraction',
'sheetFraction':'Secondary sheet fraction',
'secondaryHelixCount':'Count of secondary helix',
'secondarySheetCount':'Count of secondary sheet',
'secondaryTurnsCount':'Count of secondary turns',
'secondaryCoilCount':'Count of secondary coil',
'secondaryHelixPercentage':'Percentage of secondary helix',
'secondarySheetPercentage':'Percentage of secondary sheet',
'secondaryTurnsPercentage':'Percentage of secondary turns',
'secondaryCoilPercentage':'Percentage of secondary coil',
'cMannosylation':'C-mannosylation sites',
'proteasomalCleavage':'Proteasomal cleavages (MHC ligands)',
'nLinkedGlycosylation':'N-linked glycosylation sites',
'genericPhosphorylationSerine':'Generic phosphorylation sites of serine',
'genericPhosphorylationThreonine':'Generic phosphorylation sites of threonine',
'genericPhosphorylationTyrosine':'Generic phosphorylation sites of tyrosine',
'arginineLysinePropeptideCleavage':'Arginine and lysine propeptide cleavage sites',
'bindingRegionsDisordered':'Binding Regions in Disordered Proteins',
'mTP':'Mitochondrial targeting peptide (mTP)',
'sP':'Secretory pathway signal peptide (SP)',
'otherLocation':'Other subcellular location',
'transmembraneHelices':'Count of transmembrane helices',
'signalPeptides':'Presence of signal peptides',
'linearBcellEpitopes':'Count of linear B-cell epitopes',
'classIimmunogenicity':'Class I Immunogenicity score'}



#---------------------EUKARYOTIC BACKGROUND TASK------------------
@task
def run_eukaryotic(filename, email, features_list):    
    records=SeqIO.parse(filename, "fasta")
    featuresOutFile=open(filename+"_features.csv", "w")
    featuresOutFile.write('Sequence ID,')
    for feature in features_list[:-1]:
        featuresOutFile.write(featuresDic[str(feature)]+",")
    featuresOutFile.write(featuresDic[str(features_list[-1])])
    featuresOutFile.write('\n')
        
    for record in records:
        featuresOutFile.write(record.id+","+features(record.id, str(record.seq), features_list)+"\n")

    featuresOutFile.close()
    records.close()


    #--------------------EMAIL SENDING--------------------------------------

    if email!="":
        command = "echo 'Your SchistoProt Prediction Result is attached for job ID: "+filename+"\n\n\nKind regards,\n\nLutz Krause & Shihab Hasan\nComputational Medical Genomics Group, The University of Queensland Diamantina Institute' | mutt -a "+filename+"'_features.csv' -s 'SchistoProt Prediction Result' -- "+email
        subprocess.call(command, shell=(sys.platform!="Linux"))


    featureFile = open(filename+"_features.csv",'r')
    allLines=featureFile.readlines()
    fh, fd = allLines[:1], allLines[1:] 
    featureFile.close()
    os.remove(filename)
    os.remove(filename+"_features.csv")
    
    return (fh, fd)



#---------------------SECRETORY BACKGROUND TASK------------------
@task
def run_prokaryotic(filename, secretory_email):    
    parameter=""
    test_features=""
    test_para=""
    seqID_list=[]
    result_dict={}
    result_file=open(filename+"_result.txt", 'w')
    result_file.write("Sequence_ID\tPrediction\tScore\tProbability\n")

    records=SeqIO.parse(filename, "fasta")
    for record in records:
        hash_sequence=hashlib.md5(str(record.seq)).hexdigest()
        data=Secretory.objects.filter(sequence=hash_sequence)
        if len(data)==0:
            run_para=features(record.id, str(record.seq))
            parameter=parameter+run_para+"\n"
            seqID_list.append(record.id)
            test_features=test_features+run_para+"\n"
            test_para=test_para+record.id+","+run_para+"\n"
        else:
            for p in data:
                p.access=p.access+1
                p.time=datetime.now()
                p.save()
                test_features=test_features+p.features+"\n"
                test_para=test_para+record.id+","+p.features+"\n"
                result_dict[str(record.id)]=p.prediction+"\t"+str(p.decision)+"\t"+str(p.probability)
    records.close()
    
   #---------------------WORKING WITH SCIKIT-LEARN------------
    if parameter!="":
        parameters=StringIO(parameter)
        train_positive_file="secretory_positive.csv"
        train_negative_file="secretory_negative.csv"
        train_positive=np.genfromtxt(train_positive_file, delimiter=",")[1:,1:]
        train_negative=np.genfromtxt(train_negative_file, delimiter=",")[1:,1:]
        
        train_data=np.concatenate((train_positive,train_negative))
        train_label=np.concatenate((np.ones(len(train_positive)), np.zeros(len(train_negative))))
        test_data=np.genfromtxt(parameters, delimiter=",")
        
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
        
        train_data_scaled = min_max_scaler.fit(train_data).transform(train_data)

        test_data_scaled = min_max_scaler.fit(train_data).transform(test_data)
        

        #-------------------------------------------------------------
        clf = svm.SVC(kernel='rbf', C=12.0, gamma=2.0, probability=True)

        clf.fit(train_data_scaled, train_label)

        test_predicted = clf.predict(test_data_scaled)
        decisions = clf.decision_function(test_data_scaled)

        probabilities = clf.predict_proba(test_data_scaled)
   

        #-----------------
        duplicate_seq_remover(filename)
        fasta_rec=SeqIO.index(filename+"_nodups", "fasta")
        param=parameter.split("\n")
        i=0
        for pred in test_predicted:
            decision=str(decisions[i]).replace("[","").replace("]","")
            probability=probabilities[:,1][i]
            if pred==1.0:
                pred='Secretory Protein'
            if pred==0.0:
                pred='Non-Secretory Protein'
            test_data[i]
            if probability>0.5:
                result_dict[seqID_list[i]]=pred+"\t"+decision+"\t"+str(probability)
                p=Secretory(sequence=hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest(), prediction=pred, decision=decision, probability=probability, features=str(param[i]), access=0, time=datetime.now())    
                p.save()
            else:
                result_dict[seqID_list[i]]=pred+"\t"+decision+"\t"+str(1-probability)
                p=Secretory(sequence=hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest(), prediction=pred, decision=decision, probability=1-probability, features=str(param[i]), access=0, time=datetime.now())    
                p.save()    
            i=i+1
        fasta_rec.close()
        os.remove(filename+"_nodups")
    



    #--------------------PLOT--------------------------------------

    train_positive_file="secretory_positive.csv"
    train_negative_file="secretory_negative.csv"
    train_positive=np.genfromtxt(train_positive_file, delimiter=",")[1:,1:]
    train_negative=np.genfromtxt(train_negative_file, delimiter=",")[1:,1:]
    test_feat=StringIO(test_features)
    test_parameters=np.genfromtxt(test_feat, delimiter=",")
    train_data=np.concatenate((train_positive,train_negative))
    
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    
    train_positive_scaled = min_max_scaler.fit(train_data).transform(train_positive)
    train_negative_scaled = min_max_scaler.fit(train_data).transform(train_negative)
    test_parameters_scaled = min_max_scaler.fit(train_data).transform(test_parameters)

    x1=np.mean(train_positive_scaled, axis=0).tolist()
    x2=np.mean(train_negative_scaled, axis=0).tolist()
    x3={}

    records=SeqIO.parse(filename, "fasta")
    if test_parameters_scaled.ndim==1:
        for record in records:
            x3[record.id]=test_parameters_scaled.tolist()
            result_file.write(record.id+"\t"+result_dict[record.id]+"\n")
    else:
        j=0
        for record in records:
            x3[record.id]=test_parameters_scaled[j].tolist()
            result_file.write(record.id+"\t"+result_dict[record.id]+"\n")
            j=j+1
    result_file.close()
    records.close()


    #--------------------EMAIL SENDING--------------------------------------
    f1=open(filename+"_features.csv",'w')
    f1.write(idx+"\n"+test_para)
    f1.close()
    if secretory_email!="":
        command = "echo 'Your SchistoProt Prediction Result is attached for job ID: "+filename+"\n\n\nKind regards,\n\nLutz Krause & Shihab Hasan\nComputational Medical Genomics Group, The University of Queensland Diamantina Institute' | mutt -a "+filename+"'_result.txt' -a "+filename+"'_features.csv' -s 'SchistoProt Prediction Result' -- "+secretory_email
        subprocess.call(command, shell=(sys.platform!="Linux"))

    f = open(filename+"_result.txt",'r')
    lines = f.readlines()[1:]
    f.close()
    f2 = open(filename+"_features.csv",'r')
    allLines=f2.readlines()
    fh, fd = allLines[:1], allLines[1:] 
    f2.close()
    f1.close()
    os.remove(filename)
    os.remove(filename+"_result.txt")
    os.remove(filename+"_features.csv")
    
    return (lines, x1, x2, x3, fh, fd)



