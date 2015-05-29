'''
function to return six frame translation-to-protein of a DNA sequence

Bio.SeqUtils.six_frame_translations seems to be broken in 1.59
and removed from 
'''

from Bio.Seq import reverse_complement, translate, Seq
from Bio.Alphabet import IUPAC

class exoneratehit:
    def __init__(self,line,uid):
        '''
        parse SUGAR hit from a line
        '''
        #sugar: query_id start end strand subject_id start end strand score
        tok = line.strip().split()
        assert len(tok) == 10
        assert tok[0] == 'sugar:'
        
        self.uid = uid
        self.qid = tok[1]
        self.qstart = int(tok[2]) #already in pythonic positions
        self.qend = int(tok[3])
        self.qstrand = tok[4]
        self.sid = tok[5]
        self.sstart = int(tok[6])
        self.send = int(tok[7])
        self.sstrand = tok[8]
        self.score = int(tok[9])
        self.lines = []
        
    def append(self,line,uid):
        '''
        append a gff line to the result
        '''
        
        tok = line.strip().split('\t')[:8]
        if tok[2] == 'gene':
            tok.append('ID=%s_%d'%(self.qid,self.uid))
            self.lines.append(tok)
            tok = tok[:]
            tok[2] = 'mRNA'
            tok[8] = 'ID=mRNA_%s_%d;Parent=%s_%d'%(self.qid,self.uid,self.qid,self.uid)
            self.lines.append(tok)
        elif tok[2] == 'exon':
            tok.append('ID=exon_%s_%d;Parent=mRNA_%s_%d'%(self.qid,uid,self.qid,self.uid))
            self.lines.append(tok)
        elif tok[2] == 'cds':
            tok[2] = 'CDS'
            tok.append('ID=CDS_%s_%d;Parent=mRNA_%s_%d'%(self.qid,uid,self.qid,self.uid))
            self.lines.append(tok)

def generate_exonerate(fname):
    '''
    yield single hit records from an exonerate output file
    must have used options --showsugar --showtargetgff
    '''
    
    if type(fname) == str:
        #treat fname as filename
        f = open(fname)
    else:
        #treat fname as an open file handle
        f = fname
        
    state = 'sugar'
    for uid,line in enumerate(f):
        if state == 'sugar' and line.startswith('sugar:'):
            rec = exoneratehit(line,uid)
            state = 'findgff'
                
        elif (state == 'findgff' or state == 'finishgff') and line.startswith(rec.sid):
            rec.append(line,uid)
            state = 'finishgff'
            
        elif state == 'finishgff' and line.startswith('#'):
            yield rec
            state = 'sugar'
        
    #only close file if we opened it
    if type(fname) == str: f.close()

class blasthit:
    def __init__(self,line):
        '''
        parse one string into a blast hit record
        record match positions where start < end
        using pythonic position conventions
        add a new field reverse_match to record
        if match is reverse complement
        '''
        
        tok = line.strip().split('\t')

        flip = 0
        self.qid = tok[0]
        self.qstart = int(tok[6])
        self.qend = int(tok[7])
        self.qstrand = '+'
        if self.qstart > self.qend:
            self.qstart,self.qend = self.qend,self.qstart
            self.qstrand = '-'
        self.qstart -= 1 #make pythonic
        
        self.sid = tok[1]
        self.sstart = int(tok[8])
        self.send = int(tok[9])
        self.sstrand = '+'
        if self.sstart > self.send:
            self.sstart,self.send = self.send,self.sstart
            self.sstrand = '-'
        self.sstart -= 1 #make pythonic
        
        #record if match is reversed or not
        if flip == 1: self.reverse_match = True
        else:         self.reverse_match = False
        
        self.percent_identity = float(tok[2])
        self.alignment_length = int(tok[3])
        self.mismatches = int(tok[4])
        self.gaps = int(tok[5])
        self.e_value = float(tok[10])
        self.bit_score = float(tok[11])

def generate_blast(fname):
    '''
    this generator function iterates through a tabular blast output file
    using the default tabular output format created using eg blastn --outfmt 6
    yielding objects containing the blast hit information
    already formated, and using pythonic 0-based positions
    '''
    
    if type(fname) == str:
        #treat fname as filename
        f = open(fname)
    else:
        #treat fname as an open file handle
        f = fname
        
    for line in f:
        yield blasthit(line)
        
    #only close file if we opened it
    if type(fname) == str: f.close()
        
class gffrec:
    def __init__(self,tok,default_uid=None):
        '''
        simple line-based gff parser, does not store subfeature hierachies
        '''
        
        self.seqid = tok[0]
        self.source = tok[1]
        self.type = tok[2]
        self.start = int(tok[3])
        self.end = int(tok[4])
        self.start -= 1 #convert to pythonic position
        if tok[5] == '.': self.score = '.'
        else:             self.score = float(tok[5])
        self.strand = tok[6]
        self.phase = tok[7]
        
        self.attributes = {}
        for x in tok[8].split(';'):
            #print '|%s|'%x
            key,val = x.split('=')
            self.attributes[key] = val.split(',')

        #create convenience aliases
        if 'Parent' in self.attributes: self.parent = self.attributes['Parent'][0]
        if 'ID' in self.attributes: self.id = self.attributes['ID'][0]
        if 'Name' in self.attributes: self.name = self.attributes['Name'][0]
        
        if default_uid:
            if not 'ID' in self.attributes:
                self.attributes['ID'] = ['rjvbio.seq.gffrec.%d'%default_uid]
                self.id = self.attributes['ID'][0]

        assert self.start <= self.end
        assert self.strand in ['+','-','.','?']
        assert self.phase in ['0','1','2','.']

def generate_gff(fname,add_ids=True):
    '''
    generator function to return one line of a gff at a time
    if add_ids is True, add generic ids to items which lack ids
    '''

    if type(fname) == str:
        f = open(fname) #treat as a file name
    else:
        f = fname #treat as a file handle
        
    for uid,line in enumerate(f):
        line = line.strip()
        if line == '' or line.startswith('#'): continue
        tok = [x.strip(' ;') for x in line.split('\t')]
        assert len(tok) == 9
        
        if not add_ids: uid = None #do not add a default id to those without one
        yield gffrec(tok,uid)
        
    if type(fname) == str:
        f.close()

class gffdata:
    def __init__(self,fname,unique=[]):
        '''
        load all info from a named file or file handle
        unique lists which item types we expect to have unique ids in the file
        '''
        
        self.items = [] #items in file order
        self.ids = {}   #items stored by id
        self.kids = {}  #lists child items of each id, in file order
        
        for rec in generate_gff(fname):
            if rec.type in unique and rec.id in self.ids:
                raise Exception('duplicate Id found (%s) for item type %s'%(rec.id,rec.type))
            
            if not rec.id in self.ids: self.ids[rec.id] = []
            self.ids[rec.id].append(rec)
            self.items.append(rec)
            if not 'Parent' in rec.attributes: continue
            
            for pid in rec.attributes['Parent']:
                if not pid in self.kids: self.kids[pid] = []
                self.kids[pid].append(rec.id)

"""
def parse_gff_uids(fname):
    '''
    load all info from a gff
    enforces uniqueness of IDs
    '''
    
    items = [] #in file order
    ids = {}   #stored by id
    kids = {}  #lists child items of each id, in file order
    
    for rec in generate_gff(fname):
        assert rec.id not in ids
        ids[rec.id] = rec
        items.append(rec)
        if not 'Parent' in rec.attributes: continue
        
        for pid in rec.attributes['Parent']:
            if not pid in kids: kids[pid] = []
            kids[pid].append(rec.id)
        
    return gffdata(items,ids,kids)
"""

def generate_kmers(seq,size):
    '''
    generate all the kmers from the sequence
    '''
    
    for i in xrange(len(seq)-size+1): yield seq[i:i+size]
        

def six_frames(seq, genetic_code=1):
    '''
    input a DNA sequence
    pad to a whole number of codons if required
    probably only works for ambiguous alphabet sequences
    return the six possible protein translations
    
    '''
    
    rev = reverse_complement(seq)
        
    frames = {}
    
    for i in [0,1,2]:
        l = len(seq) - i
        j = i + l - l%3
        frames['%+d'%i] = translate(seq[i:j],genetic_code)
        frames['-%d'%i] = translate(rev[i:j],genetic_code)
        
    return frames
