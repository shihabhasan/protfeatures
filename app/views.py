import os, sys, subprocess
import hashlib, time
sys.path.append('/home/shihab/protfeatures/app')
from django.shortcuts import render_to_response, redirect
from django.http.response import HttpResponse
from django.http import HttpResponseRedirect
from django.template import RequestContext

from app.tasks import run_eukaryotic, run_prokaryotic

def index(request):
   return render_to_response('index.html')

def home(request):
   return render_to_response('index.html')

def manual(request):
   return render_to_response('help.html')

def contact(request):
   return render_to_response('contact.html', context_instance=RequestContext(request))

def thanks(request):
   name=request.POST['name']
   email=request.POST['email']
   message=request.POST['message']
   command = "echo 'Name: '"+name+"'\nEmail: '"+email+"'\nMessage: '"+message+" | mutt -s 'ProtFeatures Prediction Query' -- pharm.shihab@gmail.com"
   subprocess.call(command, shell=(sys.platform!="Linux"))
   return render_to_response('thanks.html', { 'name': name })

#-----------protfeatures PREDICTION-----------------------------

def eukaryotic_app(request):
   return render_to_response('eukaryotic_app.html', context_instance=RequestContext(request))
	
def eukaryotic_predict(request):
   features_list=request.POST.getlist('featuresBox')
   eukaryotic_email=request.POST['eukaryotic_email'].replace(" ","")
   if request.POST['eukaryotic_sequences'].replace(" ","")!="":
      filename="test_"+hashlib.md5(time.asctime()).hexdigest()
      file_in=open(filename, 'w')
      file_in.write(request.POST['eukaryotic_sequences'].replace(" >",">"))
      file_in.close()
   else:
      file = request.FILES['file']
      filename = file.name+"_"+hashlib.md5(time.asctime()).hexdigest()
      destination=open(filename,'wb+')
      for chunk in file.chunks():
         destination.write(chunk)
   task=run_eukaryotic.delay(filename, eukaryotic_email, features_list)
   task_id=task.id
   
   #return render_to_response('eukaryotic_progress.html', { 'features_list':features_list })
   return HttpResponseRedirect('/eukaryotic_progress/'+task_id, { 'task_id': task_id })

def eukaryotic_progress(request, task_id):
   result = run_eukaryotic.AsyncResult(task_id)
   while not  result.ready():
      return render_to_response('eukaryotic_progress.html', { 'task_id': task_id })
   return HttpResponseRedirect('/eukaryotic_results/'+task_id, { 'task_id': task_id })
   
def eukaryotic_results(request, task_id):
   result = run_eukaryotic.AsyncResult(task_id)
   result_data=result.get()
   fh, fd = result_data
   return render_to_response('eukaryotic_result.html', { 'fh':fh, 'fd':fd })




















