import sqlite3,os

def insert_row(cur,table,row):
    '''
    cur can be connection or cursor
    insert row into table
    guess data types from row elements
    use None for NULL value
    '''
    
    cmd = "insert into %s values ("%table
    
    for x in row:
        if x == None:
            cmd += 'NULL,'
        elif type(x) == int:
            cmd += '%d,'%x
        elif type(x) == float:
            cmd += '%e,'%x
        else:
            assert type(x) == str
            cmd += "'%s',"%x
            
    cmd = cmd[:-1] + ');'#remove last comma
    
    cur.execute(cmd)

def create_table(cur,table,spec):
    '''
    cur can be connection or cursor
    create table with the given fields
    text is default type
    
    example spec:
    name scaf bp:integer cm:real type
    '''
    
    if type(cur) == sqlite3.Connection:
        con = cur
        cur = con.cursor()
    else:
        con = None
    
    cur.execute("drop table if exists %s;"%table)
    
    tok = spec.strip().split()
    
    for i,x in enumerate(tok):
        y = x.strip().split(':')
        field = y[0]
        
        if len(y) == 1:
            _type = 'text'
        else:
            assert len(y) == 2
            _type = y[1]
    
        if i == 0:
            cur.execute("create table %s (%s %s);"%(table,field,_type))
        else:
            cur.execute("alter table %s add column %s %s;"%(table,field,_type))
            
    if con != None: con.commit()
