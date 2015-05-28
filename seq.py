'''
function to return six frame translation-to-protein of a DNA sequence

Bio.SeqUtils.six_frame_translations seems to be broken in 1.59
and removed from 
'''

from Bio.Seq import reverse_complement, translate, Seq
from Bio.Alphabet import IUPAC

class gffdata:
    def __init__(self,items,ids,kids):
        '''
        two dictionaries containing all the data from a gff file
        '''
        
        self.items = items
        self.ids = ids
        self.kids = kids

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
    '''

    if type(fname) == str:
        f = open(fname)
    else:
        f = fname
        
    for uid,line in enumerate(f):
        line = line.strip()
        if line == '' or line.startswith('#'): continue
        tok = [x.strip(' ;') for x in line.split('\t')]
        assert len(tok) == 9
        
        if not add_ids: uid = None #do not add a default id to those without one
        yield gffrec(tok,uid)
        
    if type(fname) == str:
        f.close()

def parse_gff(fname):
    '''
    load all info from a gff
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
