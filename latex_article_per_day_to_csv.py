import os, txt_mixin, copy, re

from IPython.core.debugger import Pdb

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
p_relpath = re.compile('\\\\mylink{(.*)}')
p_rating = re.compile('\\\\myrating{(.*)}')
comma_p = re.compile(' *, *')
p_sec = re.compile('\\\\section{(.*)}')
p_cite = re.compile('\\\\cite{(.*)}')


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
                                         

    def delete_empty_lines_at_begining_and_end(self, listin):
        while listin and (not listin[-1]):
            listin.pop(-1)
        while listin and (not listin[0]):
            listin.pop(0)
        return listin


    def parse_labels(self):
        self.find_subsub_titles()
        if 'Labels' not in self.sub_sub_titles:
            #quite this function
            self.labels = ''
            return
        self.label_lines = self.pop_sub_sub_section('Labels')        
        self.label_lines = filter(None, self.label_lines)
        assert self.label_lines[0] == '\\subsubsection{Labels}', "Problem with first line of label_lines"
        self.label_lines.pop(0)
        label_str = ';'.join(self.label_lines)
        label_str = comma_p.sub(';', label_str)
        self.labels = label_str


    def parse_relpath(self):
        self.find_subsub_titles()
        if 'Path' not in self.sub_sub_titles:
            #quite this function
            self.relpath = ''
            return
        self.path_lines = self.pop_sub_sub_section('Path')        
        self.path_lines = filter(None, self.path_lines)
        assert self.path_lines[0] == '\\subsubsection{Path}', "Problem with first line of path_lines"
        self.path_lines.pop(0)
        assert len(self.path_lines) == 1, "self.path_lines does not have exactly one line in it:\n" + str(self.path_lines)
        q_relpath = p_relpath.search(self.path_lines[0])
        self.relpath = q_relpath.group(1)


    def parse_rating(self):
        self.find_subsub_titles()
        if 'Rating' not in self.sub_sub_titles:
            #quite this function
            self.rating = 0
            return
        self.rating_lines = self.pop_sub_sub_section('Rating')        
        self.rating_lines = filter(None, self.rating_lines)
        assert self.rating_lines[0] == '\\subsubsection{Rating}', "Problem with first line of rating_lines"
        self.rating_lines.pop(0)
        assert len(self.rating_lines) == 1, "self.rating_lines does not have exactly one line in it:\n" + str(self.rating_lines)
        q_rating = p_rating.search(self.rating_lines[0])
        self.rating = int(q_rating.group(1))


    def find_previous_section(self):
        """Find the previous section title to use to get the month and
        year the article was read."""
        prev_ind = -1
        for ind in secinds:
            if ind > self.startind:
                #we found the first one too late
                break
            else:
                prev_ind = ind
        return prev_ind


    def get_section_month_and_year(self):
        secind = self.find_previous_section()
        secline = texlist[secind]
        q = p_sec.search(secline)
        month_year_str = q.group(1)
        month, year = month_year_str.split(' ',1)
        month = month.strip()
        year = year.strip()
        self.month_read = month
        self.year_read = int(year)


    def find_bibtex_key(self):
        cite_ind = self.clean_lines.find('\\cite{')
        cite_line = self.clean_lines[cite_ind]
        q = p_cite.search(cite_line)
        key = q.group(1)
        self.bibtex_key = key.strip()
        return self.bibtex_key
    

    def build_update_dict(self):
        mydict = {}
        mydict['relpath'] = self.relpath
	mydict['read'] = str(self.read)
	mydict['myrating'] = str(self.rating)
        mydict['mylabels'] = self.labels
	mydict['mynotes'] = self.notes
        mydict['month_read'] = self.month_read
        mydict['year_read'] = str(self.year_read)
        self.update_dict = mydict
        return mydict
    
        
    def __init__(self, linesin, startind):
        self.linesin = copy.copy(linesin)
        self.startind = startind
        self._clean_lines()
        self.parse_labels()
        self.parse_relpath()
        self.parse_rating()
        self.get_section_month_and_year()
        #pop first \\subsection line
        line0 = self.clean_lines.pop(0)
        assert line0.find('\\subsection{') == 0, "Problem with first line."
        self.find_bibtex_key()
        self.clean_lines = self.delete_empty_lines_at_begining_and_end(self.clean_lines)
        #clean up any remaining subsubsections:
        self.clean_lines.replaceallre('\\\\subsubsection{(.*)}','**\\1**')
        self.notes = ' '.join(self.clean_lines)
        self.read = 1
        
        
ind1 = 10
myparser1 = article_tex_parser(sslists[ind1], ssinds[ind1])
dict1 = myparser1.build_update_dict()

ind2 = -5
myparser2 = article_tex_parser(sslists[ind2], ssinds[ind2])
dict2 = myparser2.build_update_dict()

#for ind, mylist in zip(ssinds[10:11], sslists[10:11]):
#    myparser = article_tex_parser(mylist, ind)

#for ind, mylist in zip(ssinds[-5:-4], sslists[-5:-4]):
#    myparser = article_tex_parser(mylist, ind)


# To Do:
# - set read to 1
# - update the database
