import sys
import subprocess
workdirectory = [
"/alice/data/2016/LHC16k/000256692/pass2/PWGLF/LF_pp/1201_20190330-1815_child_7",
"/alice/data/2016/LHC16k/000257028/pass2/PWGLF/LF_pp/1201_20190330-1815_child_7",
"/alice/data/2016/LHC16k/000258359/pass2/PWGLF/LF_pp/1201_20190330-1815_child_7",
"/alice/data/2016/LHC16k/000258017/pass2/PWGLF/LF_pp/1201_20190330-1815_child_7",
"/alice/data/2016/LHC16k/000257630/pass2/PWGLF/LF_pp/1201_20190330-1815_child_7",
"/alice/data/2016/LHC16k/000256592/pass2/PWGLF/LF_pp/1201_20190330-1815_child_7"
]


for i in range(1,len(workdirectory)+1):
    digits = "{0:0=3d}".format(i)
    subprocess.call("alien_cp alien:" + workdirectory[i-1] + "/AnalysisResults.root  ./" + str(digits) + ".root", shell=True)
