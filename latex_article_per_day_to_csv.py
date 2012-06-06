import os, txt_mixin, copy, re

texpath = '/home/ryan/siue/Research/litreview/article_per_day/article_per_day.tex'

texfile = txt_mixin.txt_file_with_list(texpath)
texlist = texfile.list

ssinds = texlist.findall('\\subsection{')

sslists = []

nextinds = ssinds[1:] + [-3]

for curind, nextind in zip(ssinds, nextinds):
    curlist = texlist[curind:nextind]
    sslists.append(curlist)
    

secinds = texlist.findall('\\section{')

## \subsubsection{Labels}

## capstone design, mechatronics, education, TUES 2011

## \subsubsection{Rating}
## \myrating{3}

## \subsubsection{Path}
## \mylink{articles/education/Identifying_barriers_to_and_outcomes_of_interdisciplinarity_in_the_engineering_classroom.pdf}

p_sub_sub = re.compile('\\\\subsubsection{(.*)}')

class article_tex_parser(object):
    def find_subsub_titles(self):
        self.sub_sub_inds = self.clean_lines.findall('\\subsubsection{')
        self.sub_sub_titles = []
        for ind in self.sub_sub_inds:
            curline = self.clean_lines[ind]
            q = p_sub_sub.search(curline)
            title = q.group(1)
            self.sub_sub_titles.append(title)
            
        
    def _clean_lines(self):
        self.clean_lines = txt_mixin.txt_list(self.linesin)
        secinds = self.clean_lines.findall('\\section{')
        assert len(secinds) < 2, 'Found more than one \\section within a \\subsection'
        if len(secinds) == 1:
            self.clean_lines.pop(secinds[0])

        while not self.clean_lines[-1]:
            self.clean_lines.pop(-1)


    def pop_sub_sub_section(self, title):
        #find the start ind for title based on the index
        sub_sec_num = self.sub_sub_titles.index(title)
        start_ind = self.sub_sub_inds[sub_sec_num]
        N = len(self.sub_sub_inds)
        if sub_sec_num < (N-1):
            end_ind = self.sub_sub_inds[sub_sec_num+1]
        else:
            end_ind = None
        sub_sub_lines = copy.copy(self.clean_lines[start_ind:end_ind])
        del(self.clean_lines[start_ind:end_ind])
        return sub_sub_lines
                                         

    def get_labels(self):
        self.label_lines = self.pop_sub_sub_section('Labels')

        
    def __init__(self, linesin, startind):
        self.linesin = copy.copy(linesin)
        self.startind = startind
        self._clean_lines()
        self.find_subsub_titles()
        self.get_labels()


for ind, mylist in zip(ssinds[10:11], sslists[10:11]):
    myparser = article_tex_parser(mylist, ind)
