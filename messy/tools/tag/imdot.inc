
// = imdot.inc =======================================================
// mz_sg_20200715 (gromov@tonnerre.mpic.de)
//
// include file for imtag for production of .dot file for graphviz
// requirements: free pascal compiler (fpc) ver. 3 or higher
//
// [S.Gromov, MPIC, 2007-2020]
// ===================================================================

// - some code is borrowed from imtag v.<2 ---------------------------

// non-intermediate species list
var  nisl : array[1..max_tsl] of string;
    _nisl : integer;

// make a list of non-intermediate species
procedure make_nisl;
var s, i, j : integer;
    f : boolean;
begin
//writeln('non-intermediate species:');
_nisl:=0;
for s:=1 to _tsl do
    begin
    f:=false;
    for i:=1 to _eqs do
        with eqs[i] do
             for j:=1 to _prod do
                 if (prod[j]=tsl[s].spec) then
                    f:=true;
    if not(f) then
       begin
       inc(_nisl);
       nisl[_nisl]:=tsl[s].spec;
       end;
    end;

//writeln('#',_nisl,' species: ');
//for i:=1 to _nisl do
//    write(nisl[i],' ');
//writeln;
//writebreak;
end;

// returns true if a species is not an intermediate
function in_nisl(spec : string) : boolean;
var i : integer;
begin
in_nisl:=false;
for i:=1 to _nisl do
    if (spec=nisl[i]) then
       begin in_nisl:=true; break; end;
end;


// tagging record for a set of species
type ttagset = array[1..max_tsl] of record
     spec : nstr;                             // species name
    _prod : integer;                          // & # of reactions
     prod : array[1..max_eqs] of integer;    // reaction #
     prec : array[1..max_eqs] of string;     // who produced
     stoi : array[1..max_eqs] of real;       // & corr. stoichiometric coeff

    _sink : integer;
     sink : array[1..max_eqs] of integer;
     is2n : array[1..max_eqs] of boolean;    // sinks to nothing?
     isd2 : array[1..max_eqs] of boolean;    // sinks doubly?

    _akie : integer;                          // no of KIEs associated
     akie : array[1..max_eqs] of integer;    // +/- of abs kie r-n index in imcom
     skie : array[1..max_eqs] of integer;    // no of the isotopologue

     ifix : boolean;                          // is fixed?
     end;

// investigate the tagging information of the current mechanism and allowed tagged species
procedure tagit_v1(var d : ttagset);

var i, j : integer;
    e, r, s, p  : integer;
    qae : array[1..2] of integer;
    af : real;

begin

// make non-intermediate species list
make_nisl;

//writeln('taginfo:');

// initializing
fillchar(d,sizeof(d),0);

for s:=1 to _utsl do         // cycle through tagged species found in mech.
    with d[s] do

    begin

    // initialize record for each
    spec:=tsl[utsl[s]].spec;
   _prod:=0;
    fillchar(prod,sizeof(prod),0);
    fillchar(prec,sizeof(prec),0);

   _sink:=0;
    fillchar(sink,sizeof(sink),0);

   _akie:=0;
    fillchar(akie,sizeof(akie),0);
    fillchar(skie,sizeof(skie),0);

//  ifix:=tsl[s].ifix[..];
//  if (ifix) then continue;             // skipping the fixed species?

    for r:=1 to _eqs do // cycle through reactions list
        if (eqs[r].itag) then          // checking "dedicated for tagging" flag (additional)
           begin

           // so far skipping exchange reactions
           if (eqs[r].iiex) then
              begin
              // by convention, noting only sink for only one
              for e:=1 to 2 do
                 if (eqs[r].educ[e]=spec) then
                    begin
                    inc(_prod);
                    prod[_prod]:=r; stoi[_prod]:=1;
                    prec[_prod]:=eqs[r].educ[3-e];
                    end;
              continue; // skipping the rest
              end;

           // inspecting sink
           if ((eqs[r].educ[1]=spec) or (eqs[r].educ[2]=spec)) then
              begin
              inc(_sink);
              sink[_sink]:=r;

              // checking if sinks doubly
              isd2[_sink]:=(eqs[r].educ[1]=eqs[r].educ[2]);

              // checking if it sinks to nothing
              is2n[_sink]:=true;
              for i:=1 to eqs[r]._prod do
                  if (in_tsl(eqs[r].prod[i])) then
                     begin
                     is2n[_sink]:=false;
                     break;
                     end;
              end;

           // inspecting production
           for p:=1 to eqs[r]._prod do
               if (eqs[r].prod[p]=spec) then  // if a species is found to be produced in this r-n
                  if (eqs[r].isrc) then       // source-specified reaction, treating separately
                     for e:=1 to src[eqs[r].nsrc]._trans do     // searching through all specified transfer records
                         for j:=1 to src[eqs[r].nsrc].trans[e]._dst do     // searching in current transfer destinations list
                             if (src[eqs[r].nsrc].trans[e].dst[j].spec=spec) then
                                begin
                                inc(_prod);
                                prod[_prod]:=r;
                                stoi[_prod]:=src[eqs[r].nsrc].trans[e].dst[j].stoi*   // branching coeff. from src
                                             eqs[r].stoi[p];                          // product's self stoi
                                prec[_prod]:=src[eqs[r].nsrc].trans[e].src;
                                end
                             else
                  else  {if isrc}             // regular r-n, mass/composition weighted transfer probability
                      for e:=1 to 2 do // cycle through educts to determine whether one is tagged (or both!)
                          if (in_tsl(eqs[r].educ[e])) then
                             begin

                             // q-ty of atoms in educts
                             qae[e]:=tsl[no_tsl(eqs[r].educ[e])].qatm;
                             if (in_tsl(eqs[r].educ[3-e])) then
                                qae[3-e]:=tsl[no_tsl(eqs[r].educ[3-e])].qatm
                             else
                                 qae[3-e]:=0;

                             // accounting mass/comp. transfer balance as mass-weighted share of the products
                             af:=qae[e]/(qae[1]+qae[2]);

                             // if compound reacts with itself,
                             if (eqs[r].educ[1]=eqs[r].educ[2]) then
                                if (e=2) then  // avoiding double accounting processing 2nd educt
                                   continue
                                else           // transfer should be straight since the same educt
                                    af:=1.0;

                             // adding new product record
                             inc(_prod);
                             prod[_prod]:=r;
                             stoi[_prod]:=af*eqs[r].stoi[p];
                             prec[_prod]:=eqs[r].educ[e];
                             end;

           end;

//  writeln('-',spec,'-------------------------------------------------------------- ');

//  write('S (',_sink,'): ');
    for i:=1 to _sink do
        begin
        for j:=1 to _kie do
            if (kie[j].imec) then
               if (eqs[sink[i]].abbr=kie[j].abbr) then
                  for p:=1 to _isos do
                      if (kie[j].isot=tsl[no_tsl(spec)].isos[p]) then
                         begin
                         inc(_akie);
                         akie[_akie]:=-j;
                         skie[_akie]:=p;
//                       write('KIE[',kie[j].isot,']~');
                         end;
//      write(eqs[sink[i]].abbr,' ');
        end;
//  writeln;

//  write('P (',_prod,'): ');
    for i:=1 to _prod do
        begin

//      if (eqs[prod[i]].isrc) then
//         write('SUB:');

        for j:=1 to _kie do
            if (kie[j].imec) then
               if (eqs[prod[i]].abbr=kie[j].abbr) then
                  for p:=1 to _isos do
                      if (kie[j].isot=tsl[no_tsl(prec[i])].isos[p]) then
                         begin
                         inc(_akie);
                         akie[_akie]:=j;
                         skie[_akie]:=p;
//                       write('KIE[',kie[j].isot,']~');
                         end;

//      if not(stoi[i]=1.0) then
//         write(floattostr(stoi[i]),'*');

//      write(eqs[prod[i]].abbr);

//      write('[',prec[i],']  ');
        end;

//  writeln;

    end;

//writeln;
//writeln('tagit(',paramstr(1),',',paramstr(2),'): done');

end;

// - original imdot routines -----------------------------------------

function rem_specchar(s : ansistring) : ansistring;
begin
// avoiding control characters (like \n, \t, etc) by replacing \ to /
while (pos('\',s)>0) do
      begin
      insert('/',s,pos('\',s));
      delete(s,pos('\',s),1);
      end;
rem_specchar:=s;
end;

function rem_spaces(s : ansistring) : ansistring;
begin
// removing speces from string
while (pos(' ',s)>0) do
      delete(s,pos(' ',s),1);
rem_spaces:=s;
end;

// colors & styles for each reaction
var eqs_col : array[1..max_eqs] of string;
    eqs_sty : array[1..max_eqs] of string;

// species ranks used for species-plot
var rank_spec : array[1..max_tsl] of integer; // species rank
var prod_sink : array[1..max_tsl,1..2] of integer; // # of species producing / sink to
const rank_maxq = 50;                   // max reactions in the column

const rcm : array[0..10] of string =           // reactions colomap
    ('white',             // default (or unknown)
     'firebrick',         // oxygen, G1-G2
     'dimgray',           // nitrates, G3
     'dodgerblue3',       // 1C carbons, G41
     'khaki4',     // higher carbons, rest of G4
     'green4',            // chlorine, G6
     'lightsalmon4',      // bromine, G7
     'blueviolet',        // iod, G8
     'yellow3',           // sulphur, G9
     'darkseagreen4'{'thistle4'},           // HETs           exch nit & HET
     'orange2');          // photolysis, J


const rcm2_size = 12;
      rcm2 : array[1..rcm2_size,1..4] of string =
      (('[GJ][12][0-9]{3}[a-z]*', 'firebrick',      'react. O+H',          '1xxx+2xxx'),
       ('[GJ]3.+',                'orange2',        'react. N+H',          '3xxx'),
       ('[GJ]41.+',               'slategray3',     '1C carbons',          '41xx'),
       ('[GJ]4[2-9A-F].+',        'khaki4',         'higher HCs',          '42xx-49xx'),
       ('[GJ]6.+',                'chartreuse3',    'chlorine',            '6xxx'),
       ('[GJ]7.+',                'lightsalmon4',   'bromine',             '7xxx'),
       ('[GJ]8.+',                'blueviolet',     'iodine',              '8xxx'),
       ('[GJ]9.+',                'yellow3',        'sulphur',             '9xxx'),
       ('[GJ]10[0-9]{3,}.*',      'dimgray',        'mercury',             '10xxx'),
       ('HET[0-9].+',             'darkseagreen4',  'heterogeneous',       'HET'),
       ('IEX.+',                  'white',          'isotope exchange',    'IEX'),
       ('.+',                     'grey24',         'unknown',             '?') );

procedure make_eqs_col;
var r, i : integer;
begin

  // v0
  //  eqs_col[r]:=floattostrf(0.5-random*0.5,fffixed,7,3)+','+
  //              floattostrf(1-random*0.3,fffixed,7,3)+','+
  //              floattostrf(1-random*0.5,fffixed,7,3);
  //    eqs_col[r]:=floattostrf(1.0-0.5*(r/_eqs),fffixed,7,3)+' '+
  //              floattostrf(1.0-(r/_eqs),fffixed,7,3)+' '+
  //              floattostrf(0.5+0.5*(r/_eqs),fffixed,7,3);
  //    eqs_col[r]:=floattostrf((0.438+(r/_eqs)*(0.75-0.438)),fffixed,7,3)+' '+
  //                floattostrf((0.125+(r/_eqs)*(0.75-0.125)),fffixed,7,3)+' '+
  //                floattostrf((0.711*(1-(r/_eqs))),fffixed,7,3);
  //    eqs_col[r]:=floattostrf( r/_eqs ,fffixed,7,3)+' '+       // hue
  //                floattostrf( 0.75    ,fffixed,7,3)+' '+       // sat
  //                floattostrf( 0.75    ,fffixed,7,3);   // lig
  // v1
  //  eqs_col[r]:=rcm[0];
  //
  //    if (pos('J',eqs[r].abbr)>0) then     // photolysis
  //       eqs_col[r]:=rcm[10]
  //    else
  //    if (pos('HET',eqs[r].abbr)>0) then   // HETs
  //       eqs_col[r]:=rcm[9]
  //    else
  //    if (pos('G',eqs[r].abbr)>0) then    // gas-phase...
  //       begin
  //
  //       if ((pos('G1',eqs[r].abbr)>0) or
  //      (pos('G2',eqs[r].abbr)>0)) then eqs_col[r]:=rcm[1];  // oxygen
  //       if (pos('G3',eqs[r].abbr)>0) then eqs_col[r]:=rcm[2];   // nitrates
  //       if (pos('G41',eqs[r].abbr)>0) then eqs_col[r]:=rcm[3]   // 1C carbons
  //       else
  //       if (pos('G4',eqs[r].abbr)>0) then eqs_col[r]:=rcm[4]; // higher carbons
  //       if (pos('G6',eqs[r].abbr)>0) then eqs_col[r]:=rcm[5];     // chlorine, green
  //       if (pos('G7',eqs[r].abbr)>0) then eqs_col[r]:=rcm[6];     // bromine, brown
  //       if (pos('G8',eqs[r].abbr)>0) then eqs_col[r]:=rcm[7];     // iod, violet
  //       if (pos('G9',eqs[r].abbr)>0) then eqs_col[r]:=rcm[8];     // sulphur, yellow
  //
  //       end
  //    end;

for r:=1 to _eqs do
    begin

    eqs_col[r]:=rcm2[rcm2_size,2];
    eqs_sty[r]:='solid';

    // setting different line style to photolytic r-n
    if (upcase(eqs[r].abbr[1])='J') then eqs_sty[r]:='dashed';
    if (upcase(eqs[r].abbr[1..3])='IEX') then eqs_sty[r]:='solid, arrowhead=lvee';

    // checking patterns
    for i:=1 to rcm2_size do
        if pattmatch(eqs[r].abbr,rcm2[i,1],true) then
           begin eqs_col[r]:=rcm2[i,2]; break; end;


    end;

end;

type tscm = array[1..6,1..3] of string;     // species colormap

           // line       // fill            // font
const scm1 : tscm =
        (('black',      'lightgoldenrod',  'black'),    // from
         ('black',      'rosybrown1',      'black'),    // to
         ('black',      'darkolivegreen3', 'black'),    // bcs
         ('black',      'lightsteelblue2', 'black'),    // interm.
         ('black',      'red1',            'white'),    // def.node
         ('white', '50', '30'));  // bg, node fsize, edge fsize
      scm2 : tscm =
        (('orange',     'black',           'orange'),
         ('firebrick1', 'black',           'firebrick1'),
         ('green3',     'black',           'green3'),
         ('ivory2',     'black',           'ivory2'),
         ('gray',       'red',             'white' ),   // def.node
         ('black', '50', '40'));  // bg, node fsize, edge fsize
      scm3 : tscm =
        (('orange4'{black},      'orange',          'black'),
         ('darkolivegreen'{black}, 'olivedrab'{firebrick1},'white'{black}),
         ('black',      'lightgoldenrod4', 'black'),     // source
         ('black',      'lemonchiffon3',   'black'),     // interm.
         ('grey24',     'red',             'grey24' ),   // def.node
         ('black', '50', '40'));  // bg, node fsize, edge fsize


function spc2html_format(spec : nstr) : ansistring;
var o, s, a : ansistring;
    i, n : integer;
    is_prev_dig : boolean;
begin
o:=''; n:=no_spc(spec);
//writeln('spc2html_format: got "',spec,'" n = ',n,' _spc = ',_spc,' spc[1]="',spc[1].spec,'"');
if (n>0) then
   begin
   s:=spc[n].capt[0];
   i:=1;
   while (i<=length(s)) do
       case s[i] of
            '_', '^' : begin
                       if (s[i]='_') then a:='SUB' else a:='SUP';
                       o+='<'+a+'>';
                       if s[i+1]<>'<' then
                          begin o+=s[i+1]; i+=2; end
                       else
                           begin
                           inc(i,2);
                           while (s[i]<>'>') and (i<=length(s)) do
                                 begin o+=s[i]; inc(i); end;
                           inc(i);
                           end;
                       o+='</'+a+'>';
                       end;
            '<' : begin o+='&gt;'; i+=1; end;
            '>' : begin o+='&lt;'; i+=1; end;
            '&' : begin o+='&amp;'; i+=1; end;
            else  begin o+=s[i]; i+=1; end;
       end; {case}
   end
else
    if (spec='hv') then   // italics for hv
       o:='<I>'+spec+'</I>'
    else
        o:=spec;
spc2html_format:=o;
end;



procedure imdot_produce_dot_files(fname : string; _from, _to : string);

var d : ttagset;
    p, s, r, cp, cs, rank : integer;
    a, b, leg, prl, sign, kiestar, otheduc : ansistring;
    scm : tscm;                                      // color scheme
    rels : array[1..max_tsl,1..max_tsl+1] of string; // relations for common-relations plot
    is_dest : boolean;                               // "destruction-only" flag
    dotmech, dotspec, dotcomm, tmp : text;           // .dot files instances + template


var srcspc_list : ansistring = '';  // list of "source-only" species

function make_spec_decl(spec, specid : string) : ansistring;
var s : ansistring;
begin
s:='  '+specid+' [ label=<'+spc2html_format(spec)+'>, ';
if (pos(spec+' ',_from)>0) then
  s:=s+' shape=diamond, color='+scm[1,1]+', fillcolor='+scm[1,2]+' fontcolor='+scm[1,3]
else
  if (pos(spec+' ',_to)>0) then
    s:=s+' shape=invtriangle, color='+scm[2,1]+', fillcolor='+scm[2,2]+', fontcolor='+scm[2,3]
  else
    if (in_nisl(spec)) then
      if (is_usedspec(spec)) then
        s:=s+' color='+scm[3,1]+', fillcolor='+scm[3,2]+', fontcolor='+scm[3,3]
      else
        s:=s+' color='+scm[3,1]+', fillcolor='+scm[3,3]+', fontcolor='+scm[3,2]
    else
      if (in_tsl(spec)) then
        begin
        s:=s+' color='+scm[4,1]+', fillcolor='+scm[4,2]+', fontcolor='+scm[4,3];
        srcspc_list+=' '+spec;
        end;
s:=s+' ];';
make_spec_decl:=s;
end;


begin

_from:=trim(_from); _to:=trim(_to);

write('produce_dot_file(',fname,',',_from,' > ',_to,'): ');
_from+=' '; _to+=' ';

tagit_v1(d);

imcom_check_files_exist(['imtag.dot']);

randomize;

make_eqs_col; // prepare colormap for equations

// select colorscheme
scm:=scm3;   // cs1-3

assign(dotmech,fname+'.dotmech.dot'); rewrite(dotmech);
assign(dotspec,fname+'.dotspec.dot'); rewrite(dotspec);
assign(dotcomm,fname+'.dotcomm.dot'); rewrite(dotcomm);

a:='digraph '+tagname+' {'+_LF;
a:=a+_LF;

// copying a header from template
assign(tmp,'imtag.dot'); reset(tmp);
while not(eof(tmp)) do
      begin
      readln(tmp,prl);
      a:=a+prl+_LF;
      end;
close(tmp);

p:=0;
for s:=1 to _eqnfile do
    if (pos('wanted',eqnfile[s].line)>0) then
       begin p:=s; writeln('@$@$'); break; end;

//a+='  label = "IMTAG // CFG: '+tagname+' - ATOM: '+isoelem+'\n"+'+_LF+
//   '    "TSL: '+rem_specchar(paramstr(2))+' ['+inttostr(_utsl)+'/'+inttostr(_tsl)+'] '+
//   'EQN: '+rem_specchar(paramstr(1))+' ['+inttostr(nooftagreac)+'/'+inttostr(_eqs)+']\n"+'+_LF;
a+='  label = "IMTAG // CFG: '+tagname+' - ATOM: '+isoelem+
   ' - TAG: ['+inttostr(_utsl)+'('+inttostr(_tsl)+')/'+inttostr(_spc)+' SPC]'+
   ' ['+inttostr(nooftagreac)+'/'+inttostr(_eqs)+' EQN]\n"+'+_LF;

// finding mech info line in eqnfile
b:=GetEnvironmentVariable('xmecca_wanted4imtag');
if (b<>'') then
   a+='    "MECH: '+wraptext(rem_specchar(b),'\n',['=',' ',',','(',')','!'], 70)+'"+'+_LF;

a+='    "='+formatdatetime('DD mmm YYYY hh:nn',now)+'= [Gromov, MPI-C]";'+_LF;
a+=_LF;
a+='  // color scheme is selected in imdot.inc'+_LF;
a+='  bgcolor = '+scm[6,1]+';'+_LF;
a+=_LF;
a+='  node [color='+scm[5,1]+', fillcolor='+scm[5,2]+', fontcolor='+scm[5,3]+', style="bold, filled" ];'+_LF;
a+='  edge [color='+scm[5,1]+', fontcolor='+scm[5,3]+'];'+_LF;
a+=_LF;

writeln(dotcomm,a);
writeln(dotmech,a);
writeln(dotspec,a);

// legend (based on MECCA selection) ------------------------------------------
leg:='  // legend'+_LF;
leg+='  empty [color='+scm[6,1]+', fillcolor='+scm[6,1]+', fontcolor='+scm[6,1]+', label=" " ];'+_LF;
leg+='  intermediate [color='+scm[4,1]+', fillcolor='+scm[4,2]+', fontcolor='+scm[4,3]+'];'+_LF;
leg+='  source [color='+scm[3,1]+', fillcolor='+scm[3,2]+', fontcolor='+scm[3,3]+'];'+_LF;
//leg+='  AVEx [label = "AVE"];'+_LF;
//leg+='  AVEx -> AVEx [label = "pool average", fontcolor=gray ];'+_LF;
leg+='  LOSSx [label = "LOSS"];'+_LF;
leg+='  LOSSx -> LOSSx [label = "no educt", fontcolor=gray ];'+_LF+_LF;
leg+='  source -> intermediate;';

// v1
// leg+='  G12xx -> G12xx [label="(reac.O,H)", color='+rcm[1]+', fontcolor='+rcm[1]+'];'+_LF;
// leg+='  G3xxx -> G3xxx [label="(reac.N,H)", color='+rcm[2]+', fontcolor='+rcm[2]+'];'+_LF;
// leg+='  G41xx -> G41xx [label="(1C carbons)", color='+rcm[3]+', fontcolor='+rcm[3]+'];'+_LF;
// leg+='  G4xxx -> G4xxx [label="(higher HC)", color='+rcm[4]+', fontcolor='+rcm[4]+'];'+_LF;
// leg+='  G6xxx -> G6xxx [label="(chlorine)", color='+rcm[5]+', fontcolor='+rcm[5]+'];'+_LF;
// leg+='  G7xxx -> G7xxx [label="(bromine)", color='+rcm[6]+', fontcolor='+rcm[6]+'];'+_LF;
// leg+='  G8xxx -> G8xxx [label="(iodine)", color='+rcm[7]+', fontcolor='+rcm[7]+'];'+_LF;
// leg+='  G9xxx -> G9xxx [label="(sulphur)", color='+rcm[8]+', fontcolor='+rcm[8]+'];'+_LF;
// leg+='  HETxxx -> HETxxx [label="(heterogen.)", color='+rcm[9]+', fontcolor='+rcm[9]+'];'+_LF;
// leg+='  Jxxxx -> Jxxxx [label="Photolysis", color='+rcm[10]+', fontcolor='+rcm[10]+'];'+_LF;
//
// leg+='  G12xx [label="G1xxx/2xxx", color='+rcm[1]+', fontcolor='+rcm[1]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  G3xxx [color='+rcm[2]+', fontcolor='+rcm[2]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  G41xx [color='+rcm[3]+', fontcolor='+rcm[3]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  G4xxx [color='+rcm[4]+', fontcolor='+rcm[4]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  G6xxx [color='+rcm[5]+', fontcolor='+rcm[5]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  G7xxx [color='+rcm[6]+', fontcolor='+rcm[6]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  G8xxx [color='+rcm[7]+', fontcolor='+rcm[7]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  G9xxx [color='+rcm[8]+', fontcolor='+rcm[8]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  HETxxx [color='+rcm[9]+', fontcolor='+rcm[9]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
// leg+='  Jxxxx [color='+rcm[10]+', fontcolor='+rcm[10]+', fillcolor='+scm[6,1]+', style="bold" ];'+_LF;
//
// leg+='  { rank=max; G12xx G3xxx G41xx G4xxx G6xxx G7xxx G8xxx G9xxx Jxxxx HETxxx LOSSx intermediate source }'+_LF;
// leg+=_LF;

// v2
for r:=1 to rcm2_size do
  leg+='  rcm'+inttostr(r)+' [label="'+rcm2[r,4]+'", color='+rcm2[r,2]+', fontcolor='+rcm2[r,2]+', fillcolor = '+scm[6,1]+', style="bold"];'+_LF;
for r:=1 to rcm2_size do
  leg+='  rcm'+inttostr(r)+' -> rcm'+inttostr(r)+' [label="('+rcm2[r,3]+')", color='+rcm2[r,2]+', fontcolor='+rcm2[r,2]+'];'+_LF;
leg+='  { rank=min; ';
for r:=1 to rcm2_size do
  leg+='rcm'+inttostr(r)+' ';
leg+='LOSSx intermediate source empty }'+_LF;
leg+='//{ rank=max; LOSS }'+_LF;
leg+=_LF;

// now creating common plot (dotmech) -----------------------------------------

// species
srcspc_list:='';
for s:=1 to _utsl do
    with d[s] do
         begin
         writeln(dotmech,make_spec_decl(spec,spec));
         writeln(dotcomm,make_spec_decl(spec,spec));
         end;
srcspc_list:='// source-only species'+_LF+'//{ rank=same; '+srcspc_list+' };';
writeln(dotmech,srcspc_list);
writeln(dotcomm,srcspc_list);

prl:='';    // passed reactions list

// by species
for s:=1 to _utsl do
    with d[s] do
     begin

     // production reactions
     for r:=1 to _prod do
     begin
     sign:=':'+prec[r]+'>'+spec+'_'+eqs[prod[r]].abbr+'|';
     if (pos(sign,prl)=0) then
        begin

        if (is_kie(eqs[prod[r]].abbr)) then kiestar:='*' else kiestar:='';

        otheduc:='';
        if (eqs[prod[r]].educ[1]<>prec[r]) then otheduc:=' +'+spc2html_format(eqs[prod[r]].educ[1]);
        if ((eqs[prod[r]].educ[2]<>prec[r]) and
        (eqs[prod[r]].educ[2]<>'')) then otheduc:=' +'+spc2html_format(eqs[prod[r]].educ[2]);

        writeln(dotmech,'  ',prec[r],' -> ',spec,' ',
            '[',
            'label=<',eqs[prod[r]].abbr,kiestar,otheduc,
         // ' (#',inttostr(prod[r]),') ',
         // '[', {prec[r],'>',spec,}']',
            '>, ',
            'color=',eqs_col[prod[r]],', ',
            'style=',eqs_sty[prod[r]],', ',
            'fontcolor=',eqs_col[prod[r]],' ',
            '];');
        prl:=prl+sign;
        end;
     end;

     // sink reactions
     for r:=1 to _sink do
     begin
     is_dest:=true;

     if (is_kie(eqs[sink[r]].abbr)) then kiestar:='*' else kiestar:='';

     otheduc:='';
     if (eqs[sink[r]].educ[1]<>spec) then otheduc:=' +'+spc2html_format(eqs[sink[r]].educ[1]);
     if ((eqs[sink[r]].educ[2]<>spec) and
         (eqs[sink[r]].educ[2]<>'')) then otheduc:=' +'+spc2html_format(eqs[sink[r]].educ[2]);

     for p:=1 to eqs[sink[r]]._prod do
         begin

         sign:=':'+spec+'>'+eqs[sink[r]].prod[p]+'_'+eqs[sink[r]].abbr+'|';
         if (in_tsl(eqs[sink[r]].prod[p])) then
        begin
        is_dest:=false;
        if (pos(sign,prl)=0) then
           begin
           writeln(dotmech,'  ',spec,' -> ',eqs[sink[r]].prod[p],' ',
                   '[',
                   'label=<',eqs[sink[r]].abbr,kiestar,otheduc,
                // ' (#',inttostr(prod[r]),') ',
                // '[', {prec[r],'>',spec,}']',
                   '>, ',
                   'color=',eqs_col[sink[r]],', ',
                   'style=',eqs_sty[sink[r]],', ',
                   'fontcolor=',eqs_col[sink[r]],' ',
                   '];');
           prl:=prl+sign;
           end;
        end;
         end;
     if (is_dest) then
        writeln(dotmech,'  ',spec,' -> ','LOSS',' ',
            '[',
            'label=<',eqs[sink[r]].abbr,kiestar,otheduc,
         // ' (#',inttostr(prod[r]),') ',
         // '[', {prec[r],'>',spec,}']',
            '>, ','arrowhead=odot, ',
            'color=',eqs_col[sink[r]],', ',
            'style=',eqs_sty[sink[r]],', ',
            'fontcolor=',eqs_col[sink[r]],' ',
            '];');
     end;
     end;

writeln(dotmech);
writeln(dotmech,leg);


// now creating separated-species plot (dotspec) ------------------------------

fillchar(prod_sink,sizeof(prod_sink),0);

// by species
for s:=1 to _utsl do
    with d[s] do
         begin

         sign:='_'+inttostr(s);
         writeln(dotspec,make_spec_decl(spec,spec+sign));
         prl:='';     // using as passed-species list

         // production reactions
         for r:=1 to _prod do
             begin

             if (pos(('>'+prec[r]+'<'),prl)=0) then
                begin
                writeln(dotspec,make_spec_decl(prec[r],prec[r]+sign));
                prl:=prl+':'+'>'+prec[r]+'<';
                end;

             if (is_kie(eqs[prod[r]].abbr)) then kiestar:='*' else kiestar:='';

             otheduc:='';
             if (eqs[prod[r]].educ[1]<>prec[r]) then otheduc:=' +'+spc2html_format(eqs[prod[r]].educ[1]);
             if ((eqs[prod[r]].educ[2]<>prec[r]) and
                 (eqs[prod[r]].educ[2]<>'')) then otheduc:=' +'+spc2html_format(eqs[prod[r]].educ[2]);

             writeln(dotspec,'  ',prec[r]+sign,' -> ',spec+sign,' ',
                             '[',
                             'label=<',eqs[prod[r]].abbr,kiestar,otheduc,
                             //' (#',inttostr(prod[r]),') ',
                             //'[', {prec[r],'>',spec,}']',
                             '>, ',
                             'color=',eqs_col[prod[r]],', ',
                             'fontcolor=',eqs_col[prod[r]],' ',
                             '];');
             inc(prod_sink[s,1]);
             end;

         // sink reactions
         for r:=1 to _sink do
             begin
             is_dest:=true;

             if (is_kie(eqs[sink[r]].abbr)) then kiestar:='*' else kiestar:='';

             otheduc:='';
             if (eqs[sink[r]].educ[1]<>spec) then otheduc:=' +'+spc2html_format(eqs[sink[r]].educ[1]);
             if ((eqs[sink[r]].educ[2]<>spec) and
                (eqs[sink[r]].educ[2]<>'')) then otheduc:=' +'+spc2html_format(eqs[sink[r]].educ[2]);

             for p:=1 to eqs[sink[r]]._prod do
                 if (in_tsl(eqs[sink[r]].prod[p])) then
                    begin
                    is_dest:=false;
                    if (pos(('>'+eqs[sink[r]].prod[p]+'<'),prl)=0) then
                       begin
                       writeln(dotspec,make_spec_decl(eqs[sink[r]].prod[p],
                                       eqs[sink[r]].prod[p]+sign));
                                       prl:=prl+':'+'>'+eqs[sink[r]].prod[p]+'<';
                       end;
                    writeln(dotspec,'  ',spec+sign,' -> ',eqs[sink[r]].prod[p]+sign,' ',
                                    '[',
                                    'label=<',eqs[sink[r]].abbr,kiestar,otheduc,
                                    //' (#',inttostr(prod[r]),') ',
                                    //'[', {prec[r],'>',spec,}']',
                                    '>, ',
                                    'color=',eqs_col[sink[r]],', ',
                                    'fontcolor=',eqs_col[sink[r]],' ',
                                    '];');
                    inc(prod_sink[s,2]);
                    end;

                    if (is_dest) then
                       begin
                       writeln(dotspec,make_spec_decl('LOSS','LOSS'+sign));
                       writeln(dotspec,'  ',spec+sign,' -> ','LOSS'+sign,' ',
                                       '[',
                                       'label=<',eqs[sink[r]].abbr,kiestar,otheduc,
                                       //' (#',inttostr(prod[r]),') ',
                                       //'[', {prec[r],'>',spec,}']',
                                       '>, ','arrowhead=odot, ',
                                       'color=',eqs_col[sink[r]],', ',
                                       'fontcolor=',eqs_col[sink[r]],' ',
                                       '];');
                       inc(prod_sink[s,2]);
                       end;

                    end;

             writeln(dotspec);
             end;

// sorting layout
cp:=0; cs:=0; rank:=1; // production, sink and rank counters
fillchar(rank_spec,sizeof(rank_spec),0);
while (true) do
      begin

      r:=0;
      for p:=1 to _utsl do
      if (rank_spec[p]=0) then    // choosing not rank-assigned spec left
         begin r:=p; break; end;
      if (r=0) then break;        // no specs left, finishing

      inc(cp,prod_sink[r,1]); inc(cs,prod_sink[r,2]);
      rank_spec[r]:=rank;

      for p:=1 to _utsl do
      if (rank_spec[p]=0) then          // choosing not rank-assigned spec left
         if ((cs+d[p]._sink)<rank_maxq) and
        ((cp+d[p]._prod)<rank_maxq) then
        begin
        inc(cp,prod_sink[p,1]); inc(cs,prod_sink[p,2]);
        rank_spec[p]:=rank;
        end;

      if (cs>rank_maxq) or (cp>rank_maxq) then
     begin cp:=0; cs:=0; inc(rank); end;

      end;

a:=''; sign:='';
for r:=1 to rank do
    begin

    sign:=sign+'  { rank=same; P'+inttostr(r*3);
    for s:=1 to _utsl do
    if rank_spec[s]=r then sign:=sign+' '+tsl[utsl[s]].spec+'_'+inttostr(s);
    sign:=sign+' }'+_LF;

    a+=' P'+inttostr(r*3-1)+' -> P'+inttostr(r*3)+' -> P'+inttostr(r*3+1)+' ->';
    end;
setlength(a,length(a)-2);
writeln(dotspec,'  { node [ style=invis ]; edge [ style=invis ]; ',a,'; }');
writeln(dotspec,sign);

// now creating common-relations plot (dotcomm) -------------------------------

fillchar(rels,sizeof(rels),0);

// by species
for s:=1 to _utsl do
    with d[s] do
     begin

     // production reactions
     for r:=1 to _prod do
         if (pos(eqs[prod[r]].abbr+' ',rels[no_tsl(prec[r]),no_tsl(spec)])=0) then
            begin
            if (is_kie(eqs[prod[r]].abbr)) then rels[no_tsl(prec[r]),no_tsl(spec)]+='*';
            rels[no_tsl(prec[r]),no_tsl(spec)]+=eqs[prod[r]].abbr+' '; // ' \n';
            end;

     // sink reactions
     for r:=1 to _sink do
         begin
         is_dest:=true;

         if (is_kie(eqs[sink[r]].abbr)) then kiestar:='*' else kiestar:='';

         for p:=1 to eqs[sink[r]]._prod do
             if (in_tsl(eqs[sink[r]].prod[p])) then
                begin
                is_dest:=false;

                if (pos(eqs[sink[r]].abbr+' ',rels[no_tsl(spec),no_tsl(eqs[sink[r]].prod[p])])=0) then
                   rels[no_tsl(spec),no_tsl(eqs[sink[r]].prod[p])]+=kiestar+eqs[sink[r]].abbr+' '; // ' \n';
                end;

         if (is_dest) then
            if (pos(eqs[sink[r]].abbr+' ',rels[no_tsl(spec),_tsl+1])=0) then
               rels[no_tsl(spec),_tsl+1]+=kiestar+eqs[sink[r]].abbr+' '; // ' \n'
         end;

     end;

// relations are done, making dot

// general links
for s:=1 to _tsl do
    for p:=1 to _tsl do
        if (rels[s,p]<>'') then
           begin

           // wrapping the reactions list
           b:='<br/>' ; // line break expr. (was "\n" before html-syntax)
           rels[s,p]:=wraptext(trim(rels[s,p]), b, _delims, 16);
           // cutting last CR
         //setlength(rels[s,p],length(rels[s,p])-length(b));

           // info from > to
           // '>' changed to &gt;
           a:='&gt;' ; // "arrow" (was ">" before html-syntax)
           rels[s,p]:=spc2html_format(tsl[s].spec)+a+spc2html_format(tsl[p].spec)+':'+b+rels[s,p];

           // random colour, H-S-L coding
           a:=floattostrf(0.22-random*0.14,fffixed,7,3)+','+
              floattostrf(1.0-random*1.0,fffixed,7,3)+','+
              floattostrf(1.0-random*0.0,fffixed,7,3);

           writeln(dotcomm,'  ',tsl[s].spec,' -> ',tsl[p].spec,' ',
                           '[',
                           'label=<',rels[s,p],'>, ',
                           'color="',a,'", ',
                           'fontcolor="',a,'" ',
                           '];');

           end;

// loss reactions
for s:=1 to _tsl do
    if (rels[s,_tsl+1]<>'') then
        begin
        b:='<br/>' ; // line break expr. (was "\n" before html-syntax)
        rels[s,_tsl+1]:=wraptext(trim(rels[s,_tsl+1]), b, _delims, 16);
        // cutting last CR
      //setlength(rels[s,p],length(rels[s,_tsl+1])-length(b));
        writeln(dotcomm,'  ',tsl[s].spec,' -> LOSS ',
                           '[',
                           'label=<',    rels[s,_tsl+1]    ,'>, ',
                           'color=',rcm[1],', ',
                           'fontcolor=',rcm[1],' ',
                           '];');
        end;
writeln(dotcomm);
writeln(dotspec,leg);


// finishing ------------------------------------------------------------------

writeln(dotmech,'}');
writeln(dotspec,'}');
writeln(dotcomm,'}');

close(dotmech);
close(dotspec);
close(dotcomm);

writeln('done');

end;


(* is_prev_dig:=false;
o:='';
for i:=1 to length(spec) do
    begin
    if (spec[i] in ['0'..'9']) and not(is_prev_dig) then
       begin s+='<SUB>'; is_prev_dig:=true; end;
    if not(spec[i] in ['0'..'9']) and (is_prev_dig) then
       begin s+='</SUB>'; is_prev_dig:=false; end;
    case spec[i] of
         #10 : s+='<br/>';
         '>' : s+='&gt;';
         '<' : s+='&lt;';
         '&' : s+='&amp;';
         '"' : s+='&quot;';
         else s+=spec[i];
         end;
    end;
if (is_prev_dig) then s+='</SUB>';
dig2sub:=s;
end;
 *)
