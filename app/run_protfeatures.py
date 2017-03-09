import os, sys, subprocess
import hashlib, time
sys.path.append('/home/shihab/protfeatures/app')
from django.shortcuts import render_to_response, redirect
from django.http.response import HttpResponse
from django.http import HttpResponseRedirect
from django.template import RequestContext

from tasks import run_eukaryotic


#-----------protfeatures PREDICTION-----------------------------

	
def eukaryotic_predict(filename):
   features_list=[""]
   eukaryotic_email="pharm.shihab@gmail.com"
   task=run_eukaryotic.delay(filename, eukaryotic_email, features_list)
   task_id=task.id
   result = run_eukaryotic.AsyncResult(task_id)
   result_data=result.get()
   return result_data






















