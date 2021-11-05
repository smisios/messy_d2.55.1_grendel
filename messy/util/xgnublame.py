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

### regexp for filtering line in output of make with GNU Fortran
fname = re.compile("^.*:[0-9]+:[0-9]+:$",flags=re.IGNORECASE)
messy = re.compile("messy_",flags=re.IGNORECASE)

### Warning
ws = ['ampersand',
      'argument-mismatch',
      'c-binding-type',
      'character-truncation',
      'conversion',
      'integer-division',
      'intrinsic-shadow',
      'maybe-uninitialized',
      'return-type',
      'stringop-overflow',
      'surprising',
      'tabs',
      'target-lifetime',
      'uninitialized',
      'unused-dummy-argument',
      'unused-function',
      'unused-label',
      'unused-value',
      'unused-variable',
      'zerotrip'
      ]

allmsgs = ws
qws = []
for w in allmsgs:
    qws.append(re.compile('\[-W' + w + '\]',flags=re.IGNORECASE))
n0 = len(allmsgs)
w0 = [0] * (n0 + 1)  # + 1 for unknown

# empty dictionaries for message counters
wcount = {}

f = open('gmake.log', 'r')

while True:
    line = f.readline()
    if not line: break
    #print(line, end='')
    line = line.strip()
    if fname.search(line):
        #print(line)
        words = re.split(':',line)
        ### remove relative path and get basename
        qbase = re.split('/',words[0])
        base = qbase[-1].strip()
        ### skip three lines
        for i in range(3):
            f.readline()
            if not line: break
        msg = f.readline().strip()
        if not line: break
        #print(base + ':' + msg)
        ### analyse only messy files here
        if messy.search(base):
            if base not in wcount:
                #print('adding ' + base)
                wcount[base] = w0.copy()
            found = -1
            for q in qws:
                found = found + 1
                if q.search(msg):
                    #print(found,base,line)
                    wcount[base][found]=wcount[base][found]+1
                    break
            if found > n0:
                    wcount[base][n0+1]=wcount[base][n0+1]+1


f.close()

### dictionary contains per construction only those files, which
### have at least one Warning

print('<html>')

print('<head>')
print('   <title>')
print('    The MESSy GNU Hall of Blame')
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
print('    The MESSy GNU Hall of Blame')
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
print('  <th>warning</th>')
print('</tr>')

for f, ms in sorted(wcount.items()):
    b = sum(ms)
    if b > 0:
        print('<tr bgcolor=lightblue>')
        print('  <td><b>' + f + '</b></td>')
        print('  <td></td>')
        print('  <td></td>')
        print('</tr>')
        for i in range(len(ms)):
            if (i >= 0) and (ms[i] > 0):
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
