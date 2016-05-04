from django import template
register = template.Library()

@register.filter(name='Split')
def Split(string, index):
    return string.split("\t")[index]

@register.filter(name='split')
def split(string, token):
    return string.split(token)


