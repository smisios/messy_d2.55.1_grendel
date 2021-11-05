#!/usr/bin/env python3

import re
from datetime import datetime
import subprocess
import os

try:
    mpath = os.environ['MODULESHOME']
    cmd = [mpath + "/bin/modulecmd", "python", "list"]
    mlist = subprocess.check_output(cmd, stderr=subprocess.STDOUT).strip().decode()
except:
    mlist = '*** List of loaded modules could not be detected ***'
    pass

### where are we in our repository ? ########################################
cmd = ['git', 
       'describe', '--tags', '--always', '--dirty', '--broken']
git_version = subprocess.check_output(cmd).strip().decode()

cmd = ['git', 'branch']
git_branches = subprocess.check_output(cmd).strip().decode().splitlines()
git_branch = re.sub('\*[ ]+','',
                    [x for x in git_branches if re.search('\*', x)][0])

#cmd = ['git', 'show', '--date=iso-strict', '--format="%ad"', '--name-only']
#build_date = subprocess.check_output(cmd).strip().decode()

cmd = ['git', 'rev-parse', '--verify', 'HEAD']
git_commit = subprocess.check_output(cmd).strip().decode()

vcsrev = git_version + '_' + git_branch + '_' + git_commit

#print(vcsrev)


### get configure options ###################################################
cfopt=re.compile("^\$[ ]+\./configure .*")
f = open('config.log', 'r')
for line in f:
    line = line.strip()
    if cfopt.search(line):
        #print(line)
        confi=line
        break
f.close()

### regexp for filtering lines in output of make with NAG
warn  = re.compile("^Warning",flags=re.IGNORECASE)
quest = re.compile("^Questionable",flags=re.IGNORECASE)
messy = re.compile("messy_",flags=re.IGNORECASE)

ext   = re.compile("^Extension",flags=re.IGNORECASE)
obs   = re.compile("^Obsolescent",flags=re.IGNORECASE)
delf  = re.compile("^Deleted feature used",flags=re.IGNORECASE)

### Questionable
qs = [
    'set but never referenced',
    'is not a local variable',
    'option was not used',
    'is an unconditional EXIT statement',
    'is an unconditional RETURN statement'
    'with double precision argument and no KIND= argument returns single precision result' ]

### Warning
ws = [
    'Incompatible option setting for module',
    'explicitly imported into',
    'Unused local variable',
    'Unused intrinsic',
    'is initialised but never used',
    'but its initial value is of type',
    'Unused dummy variable ',
    'dummy argument',
    'DIMENSION specification unused because every entity has an overriding array declarator',
    'Low-precision data-value assigned to high-precision data-object',
    'referenced but never set',
    'DO loop does not loop because of unconditional EXIT statement',
    'is default initialised but never used',
    'never dereferenced',
    'Unused statement function',
    'Unused dummy procedure',
    'DO loop does not loop because of unconditional STOP statement',
    'TRANSFER to LOGICAL might produce invalid value',
    'has partly undefined result',
    'Inconsistent data type',
    'Initialisation expression for',
    'Inconsistent structure for',
    'has data type'
]

### Extensions
es = [
    'CONVERT= specifier in OPEN statement',
    'Line longer than 132 characters',
    'OUBLE COMPLEX keyword',
    'Byte count on numeric data type',
    'TAB format input',
    'inconsistent with previous argument with data type',
    'constant outside'
    ]

### Obsolescent
os = [
    'is a shared DO termination label',
    'DATA statement in executable section',
    'Statement function statement',
    'Fixed source form',
    'Arithmetic IF statement',
    'Computed GOTO statement',
    'ENTRY statement',
    'ends neither with CONTINUE nor ENDDO',
    'Assumed-length CHARACTER function'
]

### Deleted feature used
df = [
    'ASSIGN statement',
    'Assigned GOTO statement',
    'edit descriptor',
    'PAUSE statement'
]

allmsgs = ws + qs + es + os + df
qws = []
for w in allmsgs:
    qws.append(re.compile(w,flags=re.IGNORECASE))
n0 = len(allmsgs)
w0 = [0] * (n0 + 1)  # + 1 for unknown

# empty dictionaries for message counters
wcount = {}

f = open('gmake.log', 'r')
for line in f:
    #print(line, end='')
    line = line.strip()
    ### words:
    ### [0]-------------------| [1]-----------| [2]---|
    ### [Warning,Questionable]: fname, line no: message
    words = re.split(':',line)
    if ( warn.search(words[0]) or quest.search(words[0]) 
         or ext.search(words[0]) or obs.search(words[0])
         or delf.search(words[0]) ):
        ### filename
        fname = re.split(',',words[1])
        #print(fname[0])
        ### remove relative path and get basename
        qbase = re.split('/',fname[0])
        base = qbase[-1].strip()
        #print(base)
        ### analyse only messy files here
        if messy.search(base):
            if base not in wcount:
                #print('adding ' + base)
                wcount[base] = w0.copy()
            found = -1
            for q in qws:
                found = found + 1
                if q.search(words[2]):
                    #print(found,base,line)
                    wcount[base][found]=wcount[base][found]+1
                    break
            if found > n0:
                    wcount[base][n0+1]=wcount[base][n0+1]+1


f.close()

### dictionary contains per construction only those files, which
### have at least one Warning or Questionable:
### However, we want to ignore the first n message types

ignn = 0

print('<html>')

print('<head>')
print('   <title>')
print('    The MESSy NAG Hall of Blame')
print('   </title>')

print('<style>')
print('table, th, td {')
print('  border: 1px solid black;')
print('  border-collapse: collapse;')
print('}')
print('th, td {')
print('  padding: 5px;')
print('}')
print('th {')
print('  text-align: left;')
print('}')
print('</style>')

print('</head>')

print('<body>')

print('   <center>')
print('   <h1>')
print('    The MESSy NAG Hall of Blame')
print('   </h1>')
print('   </center>')

print('   <center>')
print('   <h2>')
now = datetime.now()
dt_string = now.strftime("%d.%m.%Y %H:%M:%S %Z")
print(dt_string)	
print('   </h2>')
print('   </center>')

print('<pre><code><FONT size=+2>')
print(mlist)
print('</FONT></code></pre>')

print('<h4>')
print('<FONT color=red>version   :</FONT> ' + vcsrev)
print('<p></p>')
print('<FONT color=red>configured:</FONT> ' + confi)	
print('</h4>')

print('<table style="width:100%">')
print('<tr bgcolor=black; color=white>')
print('  <th>module</th>')
print('  <th>#</th>')
print('  <th>warning/questionable/extension/obsolescent/deleted</th>')
print('</tr>')

for f, ms in sorted(wcount.items()):
    b = sum(ms[ignn:])
    if b > 0:
        print('<tr bgcolor=lightblue>')
        print('  <td><b>' + f + '</b></td>')
        print('  <td></td>')
        print('  <td></td>')
        print('</tr>')
        for i in range(len(ms)):
            if (i >= ignn) and (ms[i] > 0):
                print('<tr>')
                print('  <td></td>')
                print('  <td>', ms[i], '</TD>')
                if i <= n0:
                    print('  <td>' + allmsgs[i] + '</td>')
                else:
                    print('  <td>' + '(unknown)' + '</td>')
                print('<tr>')

print('</body>')
print('</html>')
