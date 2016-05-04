import sys, subprocess

def classIimmunogenicity(fileName):
    immunogenecity="/home/shihab/tools/immunogenicity/predict_immunogenicity.py"

    command=immunogenecity+" "+fileName
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

    lines = process.stdout.readlines()
    prediction=lines[-1].strip()
    return prediction

