import os, txt_mixin

texpath = '/home/ryan/siue/Research/litreview/article_per_day/article_per_day.tex'

texfile = txt_mixin.txt_file_with_list(texpath)
texlist = texfile.list

ssinds = texlist.findall('\\subsection{')

sslists = []

nextinds = ssinds[1:] + [-3]

for curind, nextind in zip(ssinds, nextinds):
    curlist = texlist[curind:nextind]
    sslists.append(curlist)
    
