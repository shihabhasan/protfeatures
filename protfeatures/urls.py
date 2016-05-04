"""protfeatures URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.8/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Add an import:  from blog import urls as blog_urls
    2. Add a URL to urlpatterns:  url(r'^blog/', include(blog_urls))
"""
from django.conf.urls import include, url
from django.contrib import admin

from app import views

urlpatterns = [
    url(r'^admin/', include(admin.site.urls)),
    url(r'^$', views.index, name='index'),
    url(r'^home$', views.home, name='home'),
    url(r'^help$', views.manual, name='help'),
    url(r'^contact$', views.contact, name='contact'),
    url(r'^thanks$', views.thanks, name='thanks'),
    url(r'^eukaryotic_app$', views.eukaryotic_app, name='eukaryotic_app'),
    url(r'^eukaryotic_predict$', views.eukaryotic_predict, name='eukaryotic_predict'),
    url(r'^eukaryotic_progress/(?P<task_id>[A-Za-z0-9-]+)$', views.eukaryotic_progress, name='eukaryotic_progress'),
    url(r'^eukaryotic_results/(?P<task_id>[A-Za-z0-9-]+)$', views.eukaryotic_results, name='eukaryotic_results'),
]
