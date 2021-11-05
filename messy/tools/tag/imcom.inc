
// = imcom.inc =======================================================
// mz_sg_20200715 (gromov@tonnerre.mpic.de)
//
// core code for imtag/embudget tools, version 2.6
//
// requirements: free pascal compiler (fpc) ver. 3 or higher,
//   sysutils, strutils, inifiles, classes, unix, regexpr units
//
// [S.Gromov, MPIC, 2007-2020]
// ===================================================================


// - compiler options ------------------------------------------------

// {$M+}
// {$M 50331648,33554432}
// {$M 67108864,67108864}
// {$M 120108864,120108864}
{$MODE DELPHI}
{$COPERATORS ON}


// - conditional defines ---------------------------------------------

{$DEFINE USE_REGEXP}            // use regular expressions (otherwise wildcards)
{$DEFINE USE_REGEXP_STRICT}     // apply strict matching (i.e. matched expr shoud be equal the string checked)

{$DEFINE xUSE_IF_IN_RESETPT}

{$DEFINE ADD_DUMIND}            // add dummy indices for the species listed in TSL but not present
                                // in current mech (solves missing declaration problems)

// {$DEFINE TRACDEF_ADD_H2O}    // add a fake H2O tracer record to create tagging tracers, if H2O in TSL <-- set via Makefile
// {$DEFINE TRACDEF_CHEMPROP}   // enable to use the new species/tracer info facility (chemprop+process) <-- set via Makefile


// - used libraries --------------------------------------------------

uses sysutils, strutils, inifiles, classes, unix, regexpr;


// -------------------------------------------------------------------

// confirms wildcard/regexp match of patt in str (icase=true denotes ignore case)
function pattmatch(instr, patt : string; icase : boolean) : boolean;
var n, m : boolean;
{$IFDEF USE_REGEXP}
    r : tregexpr;
{$ENDIF}
begin
if (patt[1]='!') then  // negated condition
   begin n:=true; delete(patt,1,1); end
else n:=false;
{$IFDEF USE_REGEXP}
  m:=false;
  r:=tregexpr.create;
  try
    r.expression:=patt;
    r.modifierI:=icase;
    m:=r.exec(instr);
{$IFDEF USE_REGEXP_STRICT}
    // strict matching (i.e. length of first extracted match eq. length of instr)
    if (m) then m:=(comparestr(instr,r.match[0])=0);
{$ENDIF}
  except
    on ERegExpr do
       begin
       writeln(' <!> warning: pattern ''',patt,''' is not recognised as regexp');
       m:=false;
       end;
  end;
  freeandnil(r);
{$ELSE}
  m:=iswild(instr, patt, icase);
{$ENDIF}
  if (n) then m:=not(m);
  pattmatch:=m;
end;

// - types and system chars ------------------------------------------
const _undef : real = -1E-33;                 // undefined value num. flag
const _LF = #10;                              // linefeed char
const _errstr = '<error>';                    // standart error string
const _nonestr = '<none>';                    // standart none string
const _delims = [#0,#9,#27,#32,':',';',','];  // parameters delimiters allowed
// special symbols
const _etagspsy = '?';                        // "ever tagged" in eqn ( {?#} )
const _asrssy = '&';                          // added species ref. spec ( {&spec} )
const _ifixall = '*';                         // species fixed in all classes
const _ifixone = '$';                         //               in a given class
// common type for strings, etc., species names
type nstr = string[64];

// -------------------------------------------------------------------
// ->- some fundamental physics in this pure programming
//     (masses are standard atomic, i.e. stable isotope-abundance weighted)
const _elems = 12;
const elems : array[-1.._elems] of record  // update this list in case you use missing atoms here
      name : string;
      mass : real;
      end = (
           (name:'?';  mass:1         ),  // unknown non-zero quark
           (name:'n0'; mass:01.008665 ),  // flegmatic neutron
           (name:'H';  mass:01.00794  ),  // slender hydrogen
           (name:'O';  mass:15.9994   ),  // rolling oxygen
           (name:'N';  mass:14.0067   ),  // turbo nitrogen
           (name:'C';  mass:12.0107   ),  // blackhole carbon
           (name:'Cl'; mass:35.453    ),  // sicky greeny clorine
           (name:'Br'; mass:79.904    ),  // junky bromine
           (name:'I';  mass:126.90447 ),  // silly iodine
           (name:'S';  mass:32.065    ),  // devil's sulphur
           (name:'F';  mass:18.9984   ),  // teeth fluorine
           (name:'Hg'; mass:200.59    ),  // freddy mercury
           (name:'hv'; mass:0         ),  // happy-nude-photon
           (name:'ignore'; mass:0     )); // phantom mecca bozone
// -<- voila!

// -------------------------------------------------------------------
// operation mode
{$IFDEF EMB}
const _opmode = 'budget';
{$ELSE}
const _opmode = 'tag';
{$ENDIF}

// -------------------------------------------------------------------

procedure writebreak;
begin
writeln('________________________________________________________________________________');
end;


// this here to avoid using math unit
function min(o1, o2 : longint) : longint;
begin if (o2<o1) then min:=o2 else min:=o1; end;
function max(o1, o2 : longint) : longint;
begin if (o2<o1) then max:=o1 else max:=o2; end;
function maxr(o1, o2 : real) : real;
begin if (o2<o1) then maxr:=o1 else maxr:=o2; end;


// runner showing f-n
var runner_cnt, runner_chr : integer;
procedure show_runner(freq : integer; cnt : integer = 0); inline;
  const runner = '-\|/';
//const runner = '.oOo';
begin
if (cnt<>0) then runner_cnt:=cnt else inc(runner_cnt,1);
if (runner_cnt mod max(freq,1))=0 then
   begin
   inc(runner_chr); if (runner_chr>length(runner)) then runner_chr:=1;
   write(stderr, runner[runner_chr]);
   write(stderr,#8);
   end;
end;

// issue an error and bail-out
procedure imcom_error(msg : string);
begin
writeln(' <!> error: ',msg);
halt(1);
end;

// check file for existence
procedure imcom_check_files_exist(files : array of string);
var i : byte;
begin
for i:=low(files) to high(files) do
    if not(files[i]=_nonestr) then
       if not(fileexists(files[i])) then
          imcom_error('problem, file <'+files[i]+'> is not found. stop');
end;

// --- extracts text from src located between mark1 and mark2
function imcom_ext4marks(src, markL, markR : string) : string;
var i : integer;
    s : string;
begin
if (pos(markL,src)=0) or (pos(markR,src)=0) then
   imcom_ext4marks:=''
else
    begin
    i:=pos(markL,src)+length(markL);
    s:=copy(src,i,length(src)-i+1);
    imcom_ext4marks:=copy(s,1,pos(markR,s)-1);
    end;
end;

// --- removes text from src located between mark1 and mark2, including the latter
function imcom_rem4marks(src, markL, markR : string) : string;
var s : string;
begin
s:=src;
while (pos(markL,s)>0) and (pos(markR,copy(s,pos(markL,s)+length(markL)))>0) do
   delete(s,pos(markL,s),length(markL+markR)+pos(markR,copy(s,pos(markL,s)+length(markL)))-1);
imcom_rem4marks:=s;
end;

// --- cuts the end-comment starting from expr given in markC, e.g. '//'
function imcom_rem4endcomm(src, markC : string) : string;
var s : string;
begin
s:=src;
while (pos(markC,s)>0) do
   delete(s,pos(markC,s),length(s)-pos(markC,s)+1);
imcom_rem4endcomm:=s;
end;

// function substitutes _to instead of _what in _in and gives back the result
function substr(_in, _what, _to : string) : string;
var h, l : byte;
begin
if (_in='') then exit;
h:=pos(_what,_in);
while (h>0) do
      begin
      if (_in[h-1]<>'\') then         // here '\' is a canceller for control symbols!
         begin
         delete(_in,h,length(_what));
         insert(_to,_in,h);
         inc(h,length(_to));          // advancing cursor to search in the rest of the string
         end
      else
          begin
          delete(_in,h-1,1);          // removing control char and advancing the cursor
          inc(h);
          end;
      l:=pos(_what,copy(_in,h+1,length(_in)-h));    // seaching for a new instance in the rest of str
      if (l=0) then
         h:=0
      else
          h:=h+l;
      end;
substr:=_in;
end;

// --- removes ini comments (i.e. part of the string starting from ';')
procedure imcom_preprocess_ini(fname : string);
var f, o, finc : text;      // finc for INCLUDE file
    s, p : string;
    i : integer;
    im : boolean;
begin

im:=false;     // include-mode: processing file included with {$INCLUDE <>} directive

assign(f,fname); reset(f);
assign(o,'imcom.tmp'); rewrite(o);

while not(eof(f)) do
      begin

      if (im and eof(finc)) then
         begin
         close(finc);
         im:=false;
         end;

      if (im) then readln(finc,s)
              else readln(f,s);

      s:=trim(s);
      if (s='') then continue;
      if (s[1]=';') then continue;

      if ((pos('{$INCLUDE',s)>0) and not(im)) then
         begin
         p:=imcom_ext4marks(s,'<','>');
         imcom_check_files_exist([p]);
         assign(finc,p);
         reset(finc);
         im:=true;
         writeln('imcom_preprocess_ini(',fname,'): processing $INCLUDE <',p,'>');
         continue;
         end;

      // processing lines, accounting backslashed \;, cutting end after ;
      p:=s[1];
      for i:=2 to length(s) do
          if s[i]<>';' then p+=s[i]
                       else if s[i-1]<>'\' then break
                                           else begin setlength(p,length(p)-1); p+=s[i]; end;
      // adding = to the lines to create keys
      p:=trim(p);
      if (p='') then continue;
      if (pos('=',p)=0) and (pos(']',p)=0) then p+='=';
      writeln(o,trim(p));
      end;

close(o);
close(f);

end;


// - basics, consts --------------------------------------------------

// limitations: warning! there is no implemented range-check (only by compiler)!
const max_eqs  = 5000;      // max equations qty expected
      max_spc  = 2000;      // max species qty expected
      max_prod = 40;        // max qty. of products in one original equation expected
      max_tsl  = 500;       // max no. species to tag expected
      max_isos = 32;        // max no. of isotopologues expected (or no. of classes)
      max_form = 20;        // max no. of different formers (code templates)

const eps_massbal = 1e-5;   // threshod used to detect mass-balance difference
      eps_zero = 1e-28;     // threshod used to detect zero

const wrapmargin : integer = 80;  // wrapping margin

// - tagging configuration data fields -------------------------------

// names for isotopologues & their data arrays
// nomenclature: I##E, I - isotopologue, ## - mass, E - element / example: I13C

// - general configuration parameters --------------------------------

var cmodel  : string;                      // carrier model prefix (opt.)
   _isos : integer;                        // number of isotopologues
    cfgname : string;                      // configuration name
    clsname : array[1..max_isos] of nstr;  // class names (inform. only, p.e. '12', '13' or 'oz', 'iz')
    isomass : array[1..max_isos] of real;  // isotope masses
    isoelem : string;                      // atom: C, O, N; ... if set, turns on isotope model for comp. transfer
   _isoelem : integer;                     //       corresponding entry in elems[]
    itransm : integer;                     // transfer model used (0-molecular/elemental, 1-single.sub.iso, etc.) (currently selected via isoelem)
    tagname : string;                      // tagging (or doubling) name
    kieproc : string;                      // kie-calculation (include) file former
    eqnproc : string;                      // configuration (include) file former
    cparams : string;                      // configuration parameters (via cond.defs)
   _dots    : boolean;                     // flag for doing dot-files creation
    dots    : string;                      // src/dst species to highlight in source files (given as, e.g. 'CH4 ISOP > CO')
   _roundmass : integer;                   // if set, rounds elemental masses to given dig. no (solves problems with EMAC H2O coupling)

    sisyntax, trsyntax : string;           // species index, tagged reactions naming syntax (ind_@, TAG@)
    ptsyntax, ptpsyntax, ptlsyntax : string; // passive tracers +production +loss PTs syntax
    eqnfname, spcfname : string;           // names of files created for MECCA (equation, species)
    tracdef : string;                      // tracers definition (*.tex) file name

    eqnfname_conf, spcfname_conf : string;         // names of MECCA files created for configuration (equation, species)
   _form_conf : integer;                           // # of configuration formers
    form_conf : array[1..max_form,1..2] of string; // formers (templates) filenames, output filenames
   _form_intr : integer;                           // # of inter-configuration formers
    form_intr : array[1..max_form,1..2] of string; // interconf.-formers (templates) filenames, output filenames

    nooftagreac : integer;                 // number of tagged reactions (with srcs!)
{$IFDEF USE_PT}
    ptracs_intr : ansistring;   // inter-conf passive tracers list    (':abbr>')
   _ptracs_intr : integer;      //   their count
    ptracs_conf : ansistring;   // configuration passive tracers list
   _ptracs_conf : integer;      //   their count
{$ENDIF}

    // species and equation lines that may be added to original mecca input
    addspc : array[1..max_spc] of string;
   _addspc : integer;
    addeqs : array[1..max_eqs] of string;
   _addeqs : integer;

// - tagged species list management ----------------------------------

var tsl : array[0..max_tsl] of record
        spec : nstr;                         // spec name
        qatm : integer;                      // q-ty of tagged isotope atoms (e.g. 2 for O in CH3O2) in the molecule, specified in cfg/spc
        nspc : integer;                      // corresponding # in total species list (spc)
        ival : array[1..max_isos] of nstr;   // initial value (expr. in delta or minfrac, for each minor isot.)
        ifix : array[1..max_isos] of boolean; // is a fixed spec (per class)?
        iprod, iloss : boolean;              // write production/loss PTs?
        mass : array[0..max_isos] of real;   // regular+isotopologues masses (calculated using tsl[].qatm)
        isos : array[1..max_isos] of nstr;   // isotopologues names (for tagging)
        end;
   _tsl : integer;                           // no. of TSL entries
   _fixs : integer;                          // no. of "fixed" species (updated later upon spc/eqn parsing)

{$IFDEF EMB}
// - embudget-specific declarations ----------------------------------
   var tsla : array[0..max_tsl] of record       // additional to TSL array with data
          _cats : integer;                      // # of budgeting categoriess
           ncat, ccat : array[1..max_eqs] of string;  // name and condition of each cat
           uprod, uloss : array[1..max_eqs] of boolean;
           end;
{$ENDIF}

// budgeted species list
var busl : array[1..max_tsl] of record
         spec : nstr;
         iprod, iloss : boolean;             // write production/loss PTs?
         end;
   _busl : integer;

// skipped equations list
var skel : array[1..max_tsl] of nstr;
   _skel : integer;

// used tagged species no.
var utsl : array[1..max_tsl] of integer;     // references to tsl's which are used in the curr. mech.
   _utsl : integer;

// - composition transfer-specific data ------------------------------

var src : array[1..max_eqs] of record        // array for (missing) sources and non-stoch. isotope branching specification
        abbr : nstr;                               // reaction name
        eqno : integer;                            // reaction no. (in eqs[])
         def : integer;                            // default source entry no (if no products specified)
      _trans : integer;                            // no. of sources in the current reaction
       trans : array[1..max_prod] of record        // isotope transfer information:
             src : nstr;                              // name of the "source" species,
              ib : integer;                           // if non-zero, indicates isotope-branching transfer record and equals the class #
          weight : real;                              // optional weight (probably obsolete now)
            _dst : integer;
             dst : array[1..max_prod] of record       // "destinations"
                 spec : nstr;                            // name of the "destination" species
                 stoi : real;                            // "stoichiometric" coefficient of the isotopic transfer
                 end;
             end;
        end;
   _src : integer;                              // no of source specification records

type shf_mode = (shf_none, shf_next_acc, shf_next_loss, shf_prev_acc, shf_prev_loss);

var shf : array[1..max_eqs] of record     // class-shift records
        spc, eqn : nstr;                      // species & equations (wildcarded)
        mode : shf_mode;                      // mode
        end;
   _shf : integer;                        // quantity

// - isotope-specific data -------------------------------------------

// kinetic isotope effect records
var kie, kierec : array[1..max_eqs*max_isos] of record
        abbr : nstr;                         // reaction name (abbr.)
        eqno : integer;                      // reaction no. (in eqs[])
        isot : nstr;                         // isotopologue to apply to (p.e. I12CO)
        expr : string;                       // expression (like "/alpha") where alpha is an updated KIE fractionation factor
        imec : boolean;                      // is in the meccanism?
        end;
   _kie, _kierec : integer;

// reduced mass approach KIE
var rmakie : record
           include,  exclude : array[1..max_eqs] of nstr; // names of reactions to assign RMA KIE to, and
          _include, _exclude : integer;                   // names of reactions to exclude from the previous list
    end;

// species paricipating in isotope exchange
var iesl : array[1..max_tsl] of integer; // references to tsl's which are used in the curr. mech.
   _iesl : integer;

var iex : array[1..max_eqs] of record     // array for isotope exchange reactions information
        abbr : nstr;                          // reaction name
        exspec : array[1..2] of integer;      // indices of the exchanging species in TSL
        phys : string;                        // reaction rate and other (like {%TrG}) information
        imec : boolean;                       // is in the mechanism?
        end;
   _iex : integer;                        // quantity

// --- formers (templates) management -----------------------------------------

// replacements for constants, names, etc. in {%$} templates

const _imcom_reps = 26;
       imcom_reps : array[1.._imcom_reps,1..2] of string =
                ( ('{%TIMEDATE}',''),     // time, date

                  ('{%CASE}',''),         // tagging or doubling model name
                  ('{%CMODEL}',''),       // carrier model name

                  ('{%ATOM}',''),         // atom string
                  ('{%A}',''),            // atom repeated
                  ('{%CONF}',''),         // configuration
                  ('{%C}',''),            // same

                  ('{%QSPEC}',''),        // # of tagged specs
                  ('{%NSPEC}',''),        // # of tagged specs identifier
                  ('{%QCLASS}',''),       // (same as %QISO) number
                  ('{%NCLASS}',''),       // (same as %NISO) identifier

                  ('{%QISO}',''),         // # of isotopologues
                  ('{%NISO}',''),         // # of isotopologues identifier

                  ('{%QFIX}',''),         // # of "fixed" species
                  ('{%NFIX}',''),         // # of "fixed" species identifier

                  ('{%NQATOM}',''),       // # of ist. element oatoms in molecule array identifier

                  ('{%PT}',''),           // PT syntax

                  ('{%TAG}',''),          // tagging name
                  ('{%RSIND}',''),        // array of related regular specs indentifiers
                  ('{%FSIND}',''),        // array of fixed specs identifiers
                  ('{%QCFLDIR}',''),      // # of flow directions in the (C)urrent configuration
                  ('{%QIFLDIR}',''),      // # of flow directions in all configurations, (I)nter

                  ('{%DBL}',''),          // doubling name
                  ('{%UPLIND}',''),       // indices of total unacc. prod/loss in the mech

                  ('\t',''),              // TAB symbol
                  ('\n','') );            // CR+LF symbol
       imtag_QCFLDIR = 21;
       imtag_QIFLDIR = 22;

procedure imcom_update_reps;
begin

imcom_reps[01,2]:=datetimetostr(now);
imcom_reps[02,2]:='tag';
imcom_reps[03,2]:=cmodel;
imcom_reps[04,2]:=isoelem;
imcom_reps[05,2]:=isoelem;
imcom_reps[06,2]:=cfgname;
imcom_reps[07,2]:=cfgname;

imcom_reps[08,2]:=inttostr(_utsl);
imcom_reps[09,2]:=tagname+'_'+'NTSPEC';

imcom_reps[10,2]:=inttostr(_isos);
imcom_reps[11,2]:=tagname+'_'+'NTCLASS';

imcom_reps[12,2]:=inttostr(_isos);
imcom_reps[13,2]:=tagname+'_'+'NTISO';

imcom_reps[14,2]:=inttostr(_fixs);
imcom_reps[15,2]:=tagname+'_'+'NTFIX';


imcom_reps[16,2]:=tagname+'_'+'QTATOM';
imcom_reps[17,2]:=ptsyntax;

imcom_reps[18,2]:=tagname;
imcom_reps[19,2]:=tagname+'_'+'RSIND';
imcom_reps[20,2]:=tagname+'_'+'FSIND';
imcom_reps[imtag_QCFLDIR,2]:='-1';                    // so far, updated in imtag
imcom_reps[imtag_QIFLDIR,2]:='-1';                    // so far, updated in imtag

imcom_reps[23,2]:=tagname;
imcom_reps[24,2]:=tagname+'_'+'UPLIND';

imcom_reps[25,2]:=#9;
imcom_reps[26,2]:=_LF;

end;

procedure imcom_make_reps(var s : ansistring);
var i, p : word;
begin

for i:=1 to _imcom_reps do
    begin
    p:=pos(imcom_reps[i,1],s);
    while (p>0) do
          begin
          delete(s,p,length(imcom_reps[i,1]));
          insert(imcom_reps[i,2],s,p);
          p:=pos(imcom_reps[i,1],s);
          end;
    end;

end;


// --- species list management ---------------------------------------

// --- no of the species in the tagging-species-list
function no_tsl(spec : nstr) : integer;
var i : integer;
begin
for i:=1 to _tsl do
    if (spec=tsl[i].spec) then
       begin no_tsl:=i; exit; end;
no_tsl:=0;
end;

// --- check whether the species is in the tagging-species-list
function in_tsl(spec : nstr) : boolean;
begin
in_tsl:=(no_tsl(spec)>0);
end;

// --- no of the species in the used-tagging-species-list
function no_utsl(spec : nstr) : integer;
var i : integer;
begin
for i:=1 to _utsl do
    if (spec=tsl[utsl[i]].spec) then
       begin no_utsl:=i; exit; end;
no_utsl:=0;
end;

// --- check whether the species is in the tagging-species-list
function in_utsl(spec : nstr) : boolean;
begin
in_utsl:=(no_utsl(spec)>0);
end;

// --- returns no. of class (isotopologue) if given name qualifies as a TSL spec isotopologue
function no_iso(ispec : nstr) : integer;
var k : integer;
begin
for k:=1 to _isos do
    if (pos(clsname[k],ispec)=1) then
       if in_tsl(rightstr(ispec,length(ispec)-length(clsname[k]))) then
          begin no_iso:=k; exit; end;
no_iso:=0;
end;

// --- check whether given name qualifies as a TSL spec isotopologue (i.e. as tagging spec)
function is_iso(ispec : nstr) : boolean;
begin
is_iso:=(no_iso(ispec)>0);
end;

// --- check whether the species is in the budgeting-species-list
function in_busl(spec : nstr) : boolean;      // = (if (no_busl(spec)>0))
var i : integer;
begin
for i:=1 to _busl do
    if (spec=busl[i].spec) then
       begin in_busl:=true; exit; end;
in_busl:=false;
end;

// --- no of the species in the budgeting-species-list
function no_busl(spec : nstr) : integer;
var i : integer;
begin
for i:=1 to _busl do
    if (spec=busl[i].spec) then
       begin no_busl:=i; exit; end;
no_busl:=0;
end;

// --- sources list management ---------------------------------------

// --- no of the src entry for a given reaction
function no_src(reac : nstr) : integer;
var i : integer;
begin
for i:=1 to _src do
    if (reac=src[i].abbr) then
       begin no_src:=i; exit; end;
no_src:=0;
end;

// --- returns the number of the trans with a given source spec.
function no_src_trans(src_no : integer; spec : nstr) : integer;
var i : integer;
begin
for i:=1 to src[src_no]._trans do
    if (spec=src[src_no].trans[i].src) then
       begin no_src_trans:=i; exit; end;
no_src_trans:=0;
end;

// --- returns the number of the dest. species in given src->trans
function no_trans_dst(src_no, trans_no : integer; spec : nstr) : integer;
var i : integer;
begin
for i:=1 to src[src_no].trans[trans_no]._dst do
    if (spec=src[src_no].trans[trans_no].dst[i].spec) then
       begin
       no_trans_dst:=i;
       exit;
       end;
no_trans_dst:=0;
end;

// === read the tagging configuration list ===========================
procedure imcom_read_tag_config(fname : string);

// extract the next first parameter in the line (i.e. remove it)
function ext_next(var s : ansistring) : string;
var a : string;
begin
s:=trim(s); s:=s+' ';
a:=copy(s,1,pos(' ',s));
delete(s,1,length(a)-1);
a:=trim(a);
ext_next:=a;
end;

var i, j, k, l : integer;
    s, a, b : ansistring;
    r : real;
    f : boolean;

    ini : tinifile;
    info : tstringlist;

begin

writeln('reading configuration file: ',fname);  // a ext

// preprocessing ini file for correct parsing
imcom_preprocess_ini(fname);

// reading tagging info file
ini:=tinifile.create('imcom.tmp');
info:=tstringlist.create;

// ---------------------------------------------------------
// reading TAG section first

_isos:=0; _fixs:=0;
isoelem:=''; itransm:=0;
kieproc:=''; eqnproc:='';
fillchar(isomass,sizeof(isomass),0);
{$IFDEF USE_PT}
_ptracs_conf:=0; ptracs_conf:='';
{$ENDIF}

ini.readsection('TAG',info);
if (info.count<=0) then
   imcom_error('[TAG] section information is missing. stop.');

if ( (info.indexof('configuration')<0) or
     (info.indexof('classes')<0) ) then
   imcom_error('some of (configuration,classes) required fields are missing or mistyped in [TAG] section. stop.');

// configuration name
cfgname:=extractword(1,ini.readstring('TAG','configuration',_errstr),_delims);

// classes no & names
s:=ini.readstring('TAG','classes',_nonestr);
if (s=_nonestr) then
   imcom_error('no information on tagging classes, check (classes) field in the configuration. stop.');

_isos:=wordcount(s,_delims);                  // # of classes or isotopologues
if (_isos>max_isos) then
   imcom_error('no of isotopologues/tagging classes detected ('+inttostr(_isos)+
               ') is > max_isos (='+inttostr(max_isos)+'). you can increase this value in the code. stop.');

s:=ini.readstring('TAG','classes',_errstr);
for i:=1 to _isos do
    clsname[i]:=extractword(i,s,_delims);

// optional parameters
cmodel:=extractword(1,ini.readstring('TAG','cmodel',_nonestr),_delims);         // carrier model (added prefix)

// optional parameters for isotopes
isoelem:=extractword(1,ini.readstring('TAG','iatom',_nonestr),_delims);         // atom (optional)
kieproc:=extractword(1,ini.readstring('TAG','kieproc',_nonestr),_delims);       // kie-processing file (optional if no KIE specified)
eqnproc:=extractword(1,ini.readstring('TAG','eqnproc',_nonestr),_delims);       // eqn-processing file (optional if no KIE specified)
s:=ini.readstring('TAG','imasses',_nonestr);
if ((s<>_nonestr) and (s<>'')) then
   for i:=1 to _isos do
       begin
       isomass[i]:=strtofloatdef(extractword(i,s,_delims),-1.0);                // isotope masses
       if (isomass[i]=-1.0) then
          imcom_error('problem recognizing #'+inttostr(i)+' isotope mass ('+extractword(i,s,_delims)+'). stop');
       end;

// species to budget
_busl:=0; fillchar(busl,sizeof(busl),0);
s:=ini.readstring('TAG','budget',_nonestr);
if ((s<>_nonestr) and (s<>'')) then
   for i:=1 to wordcount(s,_delims) do
       begin
       inc(_busl);
       busl[_busl].spec:=extractword(i,s,_delims);
       busl[_busl].iprod:=false;
       busl[_busl].iloss:=false;
       end;

// equations to skip
_skel:=0; fillchar(skel,sizeof(skel),0);
s:=ini.readstring('TAG','eqnskip',_nonestr);
if ((s<>_nonestr) and (s<>'')) then
   for i:=1 to wordcount(s,_delims) do
       begin
       inc(_skel);
       skel[_skel]:=extractword(i,s,_delims);
       end;

// create dot files?
_dots:=false; dots:='';
s:=ini.readstring('TAG','dots',_nonestr);
if ((s<>_nonestr) and (s<>'')) then
   begin
  _dots:=true;
   dots:=s;
   end;

// round up isotopologues' masses in tracdef (by def up to 8 dig.)
s:=ini.readstring('TAG','roundmass',_nonestr);
_roundmass:=strtointdef(s,8);

// user configuration parameters (cond.defs.)
cparams:='';
info.clear;
ini.readsection('CPARAMS',info);
for i:=0 to info.count-1 do
    cparams+=info[i]+' ';

// internal naming (so far)
tagname:='tag_'+cfgname;        // name of tagging, p.e. tag_IC

// un peu d'info
if cmodel='' then cmodel:=_nonestr;
_isoelem:=-1;
if isoelem='' then
   isoelem:=_nonestr
else // searching for corresponding element for
    begin
    for i:=1 to _elems do
        if (upcase(isoelem)=upcase(elems[i].name)) then
           _isoelem:=i;
    if (_isoelem=-1) then
       writeln(' <!> warning: selected element ('+isoelem+') element is not found in the reference table (elems[]). <may be normal>');
    if (isoelem<>_nonestr) then itransm:=1 else itransm:=0;   // setting for isotope-spec. transfer, if isoelem is given (upgradable)
    end;
if kieproc='' then kieproc:=_nonestr;
if eqnproc='' then eqnproc:=_nonestr;

writeln('----- configuration: ',cfgname,' -----');
writeln('      carrier model: ',cmodel);
writeln('      tagging model: ',tagname);
  write('            classes: ',_isos,' [ ');
for i:=1 to _isos do
    write(clsname[i],' ');
writeln(']');
writeln('           eqn-proc: ',eqnproc);
writeln('         parameters: ',cparams);
  write('  species to budget: '); for i:=1 to _busl do write(busl[i].spec,' '); writeln;
  write('  equations to skip: '); if (_skel=0) then write(_nonestr) else for i:=1 to _skel do write(skel[i],' '); writeln;
  write('   create dot files: '); if _dots then write('yes (',dots,')') else write('no'); writeln;
writeln('-- isotope-specific: --' );
  write('      transfer mode: ');
  if (itransm=1) then writeln('single-substituted isotopologues')
                 else writeln('molecular/elemental (not isotopic)');
writeln('              iatom: ',isoelem,' [',_isoelem,']');
  write('            imasses: ');
if (isomass[1]>0.0) then
   begin
   write('[ ');
   for i:=1 to _isos do
      write(floattostr(isomass[i]),' ');
   writeln('] (',elems[_isoelem].mass,')');
   end
else
    writeln('none specified');
writeln('     tracdef masses: rounded up to ',_roundmass,' digit after comma');
writeln('           kie-proc: ',kieproc);
writeln;

if (cmodel=_nonestr) then cmodel:='';
if (isoelem=_nonestr) then isoelem:='';

// configuration is known by now, updating replacements for templates
imcom_update_reps;

// ---------------------------------------------------------
// MECCA-specific files

// modified MECCA equation file name
eqnfname:=extractword(1,ini.readstring('MECCA','eqn',(cmodel+'_'+tagname+'_mecca.eqn')),_delims);
// modified MECCA additional species file
spcfname:=extractword(1,ini.readstring('MECCA','spc',(cmodel+'_'+tagname+'_mecca.spc')),_delims);
// species' index syntax (opt., default is ind_@)
sisyntax:=extractword(1,ini.readstring('MECCA','sisyntax','ind_@'),_delims);
// passive tracers (PTs) naming syntax (opt., default is XPT@)
ptsyntax:=extractword(1,ini.readstring('MECCA','ptsyntax','XPT@'),_delims);
   // define prod/loss PTs particular syntax here
   ptlsyntax:=extractword(1,ini.readstring('MECCA','ptlsyntax',substr(ptsyntax,'@','L')+'@'),_delims);
   ptpsyntax:=extractword(1,ini.readstring('MECCA','ptpsyntax',substr(ptsyntax,'@','P')+'@'),_delims);
// tagged reac. naming syntax (opt., default is TAG@)
trsyntax:=extractword(1,ini.readstring('MECCA','trsyntax','TAG@'),_delims);
// tracers def. tex file
tracdef:=extractword(1,ini.readstring('MECCA','tracdef',(cmodel+'_'+tagname+'.tex')),_delims);

// configuration eqn/spc files
// equation file name
eqnfname_conf:=extractword(1,ini.readstring('MECCA','cfgeqn',(cmodel+'_'+tagname+'.eqn')),_delims);
// additional species file
spcfname_conf:=extractword(1,ini.readstring('MECCA','cfgspc',(cmodel+'_'+tagname+'.spc')),_delims);

imcom_make_reps(sisyntax);
imcom_make_reps(ptsyntax);
imcom_make_reps(eqnfname); imcom_make_reps(spcfname);
imcom_make_reps(tracdef);

imcom_make_reps(eqnfname_conf); imcom_make_reps(spcfname_conf);

// un peu d'info
writeln('MECCA-specific parameters:');
writeln('    spec.ind syntax: ',sisyntax);
writeln('         PTs syntax: ',ptsyntax);
writeln(' tagged reac syntax: ',trsyntax);
writeln('     equation files: ',eqnfname,' / ',eqnfname_conf);
writeln('      species files: ',spcfname,' / ',spcfname_conf);
writeln;

// ---------------------------------------------------------
// code former (template) files

_form_conf:=0;

// reading configuration formers
a:='CODE:CONF';
info.clear;
ini.readsection(a,info);

if (info.count>max_form) then
   imcom_error('too many code templates to process. see read_tag_config and max_form. stop.');

for i:=0 to info.count-1 do
 if (info[i]<>'') then   // necessary condition to avoid empty keys of occasional trash in cfg
    begin
    inc(_form_conf);           // number of templates

    form_conf[_form_conf,1]:=trim(info[i]);          // former filename

    // reading postfix
    s:=ini.readstring(a,info[i],_errstr);

    // making the output filename   (p.e. 'messy_mecca_tag_IC_box.f90', where postfix is '_box')
    form_conf[_form_conf,2]:=trim(extractword(1,s,_delims));

    imcom_make_reps(form_conf[_form_conf,1]);
    imcom_make_reps(form_conf[_form_conf,2]);

    end;

// un peu d'info
writeln('#',_form_conf,' configuration code template(s):');
for i:=1 to _form_conf do
    writeln(form_conf[i,1]:25,' -> ',form_conf[i,2]);


// reading inter-configuration formers
a:='CODE';
info.clear;
ini.readsection(a,info);

writeln('# of inter-configuration code template(s):');

if ((_form_intr+info.count)>max_form) then
   imcom_error('too many code templates to process. see read_tag_config and max_form. stop.');

for i:=0 to info.count-1 do
 if (info[i]<>'') then   // necessary condition to avoid empty keys of occasional trash in cfg
    begin
    inc(_form_intr);           // number of templates

    form_intr[_form_intr,1]:=trim(info[i]);          // former filename

    // reading postfix
    s:=ini.readstring(a,info[i],_errstr);

    // making the output filename   (p.e. 'messy_mecca_tag_IC_box.f90', where postfix is '_box')
    form_intr[_form_intr,2]:=trim(extractword(1,s,_delims));

    imcom_make_reps(form_intr[_form_intr,1]);
    imcom_make_reps(form_intr[_form_intr,2]);

    write(form_intr[_form_intr,1]:25,' -> ',form_intr[_form_intr,2]);

    // checking tag inter-conf formers
    for j:=1 to _form_intr-1 do
        if ((form_intr[j,1]=form_intr[_form_intr,1]) and
            (form_intr[j,2]=form_intr[_form_intr,2])) then
           begin
           dec(_form_intr);  // removing last, which is detected to be duped
           write(' (detected as duplicate)');
           break;
           end;

    writeln;
    end;

writeln;

// ---------------------------------------------------------
// reading species section (SPC)

_tsl:=0; _utsl:=0;
fillchar(tsl,sizeof(tsl),0);
fillchar(utsl,sizeof(utsl),0);

tsl[0].spec:='<not found>';

ini.readsection('SPC',info);
if (info.count<=0) then
   imcom_error('[SPC] section information is missing. stop.');

info.clear;
ini.readsection('SPC',info);
for i:=0 to info.count-1 do
 if (info[i]<>'') then   // necessary condition to avoid empty keys of occasional trash in cfg
    if (not(in_tsl(extractword(1,info[i],_delims)))) then // avoid mistyped double
       begin
       inc(_tsl);                         // number of specs. to be tagged

       tsl[_tsl].spec:=trim(extractword(1,info[i],_delims)); // species name

       // # of isot. elem. atoms in the molecule
       j:=strtointdef(extractword(2,info[i],_delims),0);
       tsl[_tsl].qatm:=j;                 // # of isot. elem. atoms in the molecule

       // reading species initial compostition data
       s:=ini.readstring('SPC',info[i],_errstr);

       // number of values in the line, filling in backward-direction
       //   i.e. for isotopes, I16 I17 I18 reading first 18th, then 17th, 16th isn't necessary
       //        for fractions, F1 F2 F3 are read 3 - 2 - 1 completing all, or 3 - 2, leaving the 1st out
       k:=min(wordcount(s,_delims),max(_isos-itransm,1)); // bugfix 151104
       for j:=_isos downto (_isos-k+1) do
           begin
           tsl[_tsl].ival[j]:=extractword(j-(_isos-k),s,_delims);
        // indicated as fixed-in-class? (spec. character)
           if (pos(_ifixone,tsl[_tsl].ival[j])>0) then
              tsl[_tsl].ifix[j]:=true;
        // remove _ifixone character(s)
           tsl[_tsl].ival[j]:=substr(tsl[_tsl].ival[j],_ifixone,'');
           end;

//     for j:=min(_isos,2) to _isos do      // initial value (delta or frac), (opt., def 0)
//         tsl[_tsl].ival[j]:=extractword(max(j-1,1),s,_delims);  // for one class read at least one value, for more read values for minors

       // is fixed species (all classes)?
       if (pos(_ifixall,s)>0) then
          for j:=1 to _isos do
              tsl[_tsl].ifix[j]:=true;
       // add prod/loss PTs?
       tsl[_tsl].iprod:=(pos('P',s)>0);
       tsl[_tsl].iloss:=(pos('L',s)>0);
       end
     else
         writeln(' <!> warning: double entry for ',info[i]);

// filling isotopologues names
// -- make corrections to names here
for i:=1 to _tsl do
    for j:=1 to _isos do
        tsl[i].isos[j]:=clsname[j]+tsl[i].spec;

// un peu d'info
writeln('#',_tsl,' species to '+_opmode+': (initialisation units to be defined in the code, ''?'' denotes that atoms count is to be picked from SPC)');
  write('name':10,' (#',isoelem,') / ');
for j:=1 to _isos do
    write(clsname[j],isoelem,' ');
writeln;

for i:=1 to _tsl do
    with tsl[i] do
        begin
        write(spec:10,' (');
        if (qatm>0) then write(qatm) else write('?');
        write(isoelem,') / ');
        f:=false;
        for j:=1 to _isos do
            begin
            if (ifix[j]) then begin write('$'); f:=true; end;
            write(ival[j],' '); 
            end;
        if (f) then write(' (some or all species indicated as fixed composition)');
        if (iprod) then write(' (adding production PT)');
        if (iloss) then write(' (adding loss PT)');
        writeln;
        end;
writeln;


// ---------------------------------------------------------
// reading KIE info section

_kie:=0; _kierec:=0;
fillchar(kie,sizeof(kie),0);

info.clear;
ini.readsection('KIE',info);
for i:=0 to info.count-1 do
 if (info[i]<>'') then   // necessary condition to avoid empty keys due to occasional trash in cfg
    begin
    inc(_kierec);
    kierec[_kierec].abbr:=extractword(1,trim(info[i]),_delims);  // reaction abbr.
    kierec[_kierec].isot:=extractword(2,trim(info[i]),_delims);  // isotopologue

    // expression = whole string
    s:=ini.readstring('KIE',info[i],_errstr);
    kierec[_kierec].expr:=trim(s);

    kierec[_kierec].eqno:=-1;               // eq. no in eqs[] is not known yet
    kierec[_kierec].imec:=false;            // by default, not present in the mech
    end;

writeln('#',_kierec,' KIE-specific record(s): ');
for i:=1 to _kierec do
    with kierec[i] do
         writeln(' ',abbr,' : for ',isot,' [ ',expr,' ]');
writeln;

// and RMA KIE information
fillchar(rmakie,sizeof(rmakie),0);

info.clear;
ini.readsection('KIE:RMA',info);

if (info.count>0) then
   begin

   // processing include-list
   a:=trim(ini.readstring('KIE:RMA','include',''));
   if (a<>'') then
      for i:=1 to wordcount(a,_delims) do
          with rmakie do
               begin
               inc(_include);
               include[_include]:=extractword(i,a,_delims);
               end;

   // and exclude-list
   a:=trim(ini.readstring('KIE:RMA','exclude',''));
   if (a<>'') then
      for i:=1 to wordcount(a,_delims) do
          with rmakie do
               begin
               inc(_exclude);
               exclude[_exclude]:=extractword(i,a,_delims);
               end;

   with rmakie do
        if (_include>0) or (_exclude>0) then
           begin
           writeln('RMA KIE info:');
             write(' #',_include,' to include ( ');
           for i:=1 to _include do
               write(include[i],' ');
           writeln(')');
             write(' #',_exclude,' to exclude ( ');
           for i:=1 to _exclude do
               write(exclude[i],' ');
           writeln(')');
           end;
   end
else
    writeln('RMA KIE: no information');

writeln;

// ---------------------------------------------------------
// reading IEX info section

_iex:=0;
fillchar(iex,sizeof(iex),0);

info.clear;
ini.readsection('IEX',info);
for i:=0 to info.count-1 do
 if (info[i]<>'') then   // necessary condition to avoid empty keys due to occasional trash in cfg
    begin
    inc(_iex);
    iex[_iex].abbr:=trim(info[i]);  // reaction abbr.

    // expression = whole string
    s:=ini.readstring('IEX',info[i],_errstr);
    for j:=1 to 2 do
        iex[_iex].exspec[j]:=no_tsl(extractword(j,s,_delims));  // exchanging species
    with iex[_iex] do
         if ((exspec[1]*exspec[2])=0) then    // in_tsl returns 0 if the spec is not found in TSL
            begin
            dec(_iex);
            writeln(' <!> warning: parsing [IEX] ',abbr:6,' for ',extractword(1,s,_delims),
                    ' <-> ',extractword(2,s,_delims),' : some of species are not found in TSL');
            end
         else
             begin
             phys:=trim(copy(s,pos(':',s)+1,length(s)-pos(':',s)));
             if (pos(';',phys)=0) then  // adding ';' to the rate expr. (limitation of tinifile unit)
                phys+=';';
             imec:=false;              // by default, not present in the mech
             end;
    end;

writeln('#',_iex,' isotope exchange reaction(s): ');
for i:=1 to _iex do
    with iex[i] do
         writeln(abbr:6,': ',tsl[exspec[1]].spec,' <-> ',tsl[exspec[2]].spec,' ',phys);
writeln;

// ---------------------------------------------------------
// reading sources specification section (SRC)

write('source specification record(s): ');

_src:=0;
fillchar(src,sizeof(src),0);

a:='';
info.clear;
ini.readsection('SRC',info);
for i:=0 to info.count-1 do
 if (info[i]<>'') then   // necessary condition to avoid empty keys of occasional trash in cfg
    begin
    // reaction abbr.
    s:=extractword(1,trim(info[i]),_delims);

    j:=no_src(s);           // if already in the list, getting the entry no
    if (j=0) then
       begin                // otherwise creating a new entry
       inc(_src);
       src[_src].abbr:=s;
       j:=_src;
       end;

    with src[j] do
         begin

         // getting source species name and opt. weight
         s:=trim(extractword(2,info[i],_delims));
         if (s[1] in ['.','-','+','0'..'9']) then     // determine optional coefficient
            begin
            r:=strtofloatdef(s,1.0);                  // storing stoi. coefficient value
            s:=trim(extractword(3,info[i],_delims));  // reading following species name
            end
         else
             r:=1.0;

         if not(in_tsl(s)) and not(is_iso(s)) then    // not a tagged species, not an IB record - skipping
            begin a+=s+' '; continue; end;

         k:=no_src_trans(j,s);
         if (k=0) then
            begin
            inc(_trans);                              // new source entry
            trans[_trans].src:=s;
            trans[_trans].weight:=r;
            k:=_trans;
            trans[_trans].ib:=no_iso(s);              // if >0 = # of class for isotope branching (IB) record
            end;

         // destination species
         s:=ini.readstring('SRC',info[i],_errstr);
         s:=trim(s);

         l:=1;
         if (s='') and (trans[_trans].ib=0) then      // indicates that all products are from current source
            def:=_trans                               // not reading further
         else
             while (l<=wordcount(s,_delims)) do
                   begin
                   b:=trim(extractword(l,s,_delims));         // b - destination spec. (isotopologue if IB) name

                   if (b[1] in ['.','-','+','0'..'9']) then   // determine optional stoi. coefficient
                      begin
                      r:=strtofloatdef(b,1.0);                // storing coefficient value
                      inc(l);
                      b:=trim(extractword(l,s,_delims));      // reading species name that follows
                      end
                   else 
                       r:=1.0;

                   if not(in_tsl(b)) then
                      a+=b+' '
                   else
                       if (no_trans_dst(j,k,b)=0) then
                          with trans[_trans] do
                               begin
                               inc(_dst);
                               dst[_dst].spec:=b;
                               dst[_dst].stoi:=r;
                               end;
                   inc(l);
                   end;
         end;

    src[_src].eqno:=0;   // equation no. in eqs[], not initialized so far
    end;

if _src>0 then writeln('#',_src) else writeln('no');

for i:=1 to _src do
  with src[i] do
    begin
    write(abbr:6,': (',_trans:2,', default: ');
    if (def>0) then write(trans[src[i].def].src) else write('-');
    write(' ) : ');
    // transfers
    for j:=1 to src[i]._trans do
      with trans[j] do
        begin
        if (ib>0) then
           write('(IB) ');
        if (weight<>1.0) then
           write(floattostr(weight)+' ');
        write(src,' ->');
        for k:=1 to _dst do
            begin
            write(' ');
            if (dst[k].stoi<>1.0) then
               write(floattostr(dst[k].stoi)+' ');
//             write(formatfloat('##0.#####',src[i].trans[j].dst[k].stoi)+' ');
            write(dst[k].spec);
            end;
        if (j<_trans) then write(', ');
        end;
    writeln;
    end;
if (a<>'') then
   writeln('following species were not recognized (not in TSL?): ',a);
writeln;

// ---------------------------------------------------------
// reading additionally specified species/equations

// species
_addspc:=0; fillchar(addspc,sizeof(addspc),0);
info.clear; ini.readsection('ADD:SPC',info);
for i:=0 to info.count-1 do
    begin
    inc(_addspc);
    s:=ini.readstring('ADD:SPC',info[i],'');
    addspc[_addspc]:=trim(info[i]);
    if (s<>'') then addspc[_addspc]+=' = '+s;
    end;
if (_addspc>0) then
   begin
   writeln('additional SPC lines specified: ');
   for i:=1 to _addspc do
       writeln(addspc[i]);
   writeln;
   end;

// equations
_addeqs:=0; fillchar(addeqs,sizeof(addeqs),0);
info.clear; ini.readsection('ADD:EQN',info);
for i:=0 to info.count-1 do
    begin
    inc(_addeqs);
    s:=ini.readstring('ADD:EQN',info[i],'');
    addeqs[_addeqs]:=trim(info[i]);
    if (s<>'') then addeqs[_addeqs]+=' = '+s;
    end;
if (_addeqs>0) then
   begin
   writeln('additional EQN lines specified: ');
   for i:=1 to _addeqs do
       writeln(addeqs[i]);
   writeln;
   end;

// ---------------------------------------------------------
// reading classes-shift information

write('classes-shift record(s): ');

_shf:=0;
fillchar(shf,sizeof(shf),0);

a:='';
info.clear;
ini.readsection('SHF',info);
for i:=0 to info.count-1 do
    if (wordcount(info[i],_delims)=2) then
       begin
       inc(_shf);
       //writeln('reading :: ',info[i]);
       with shf[_shf] do
            begin
            spc:=extractword(1,info[i],_delims);
            mode:=shf_none;
            a:=extractword(2,info[i],_delims);
            if (a='+' ) then mode:=shf_next_loss;
            if (a='++') then mode:=shf_next_acc;
            if (a='-' ) then mode:=shf_prev_loss;
            if (a='--') then mode:=shf_prev_acc;
            if (mode=shf_none) then
               begin
               writeln('error parsing record (expecting <spec> <mode> = <reac>): ',info[i]);
               dec(_shf);
               end
            else
                begin
                a:=trim(ini.readstring('SHF',info[i],''));
                k:=wordcount(a,_delims);
                for j:=1 to k do
                    begin
                    shf[_shf].eqn:=extractword(j,a,_delims);
                    if (j<k) then
                       begin
                       inc(_shf);
                       shf[_shf].spc:=shf[_shf-1].spc;
                       shf[_shf].mode:=shf[_shf-1].mode;
                       end;
                    end;
                end;
            end;
       end
    else
        writeln('error parsing record (expecting <spec> <mode> = <reac>): ',info[i]);

if (_shf>0) then
   begin
   writeln(_shf);
   if (itransm>0) then writeln(' <!> warning: selected isotope composition transfer model is likely not compatible with shifts! know what you are doing.');
   for i:=1 to _shf do
       with shf[i] do
            begin
            case mode of
                 shf_next_acc  : a:='next with accumulation';
                 shf_next_loss : a:='next with loss';
                 shf_prev_acc  : a:='previous with accumulation';
                 shf_prev_loss : a:='previous with loss';
            end; { case }
            writeln(' ',spc,' in ',eqn, ', shift mode: ',a);
            end;
   end
else
    writeln('none');

writeln;

// ---------------------------------------------------------

info.destroy;
ini.destroy;

writeln('read_tag_config: done');

writebreak;
end;


// - management des equations ----------------------------------------

// equation file by lines
var eqnfile : array[1..3*max_eqs] of record
            line : ansistring;        // line read
            iseq : boolean;           // "is an equation record?" flag
            eqno : word;              // equation # (supposed to coincide with internal # in mecca)
            end;
   _eqnfile : word;                      // # of lines in eqnfile

var spcfile : array[1..3*max_spc] of record
            line : ansistring;        // line read
            ispc : boolean;           // "is a species record?" flag
            nspc : word;              // species # (k.a. if supposed to coincide with internal # in mecca)
            end;
   _spcfile : word;                      // # of lines in spcfile

// species data-array
var spc : array[0..max_spc] of record        // zero entry is for junk
        spec : nstr;                         // name
        comp : array[0..max_isos] of nstr;   // composition
        capt : array[0..max_isos] of string; // {@caption}, after the ;
        qatm : integer;                      // q-ty of tagged isotope atoms (e.g. 2 for O in CH3O2) in the molecule, picked from spc
        mass : real;                         // calculated with elem data (see below)
        tslno : integer;                     // corresponding # in TSL
        asrs : nstr;                         // added species refspec (indicates that imtag should add a trac.def.)
        end;
   _spc : integer;                    // number of species found in spc file

// equations data-array
var eqs : array[1..max_eqs] of record
        eqnf : integer;               // no of entry in eqnfile[]
        abbr : nstr;                  // abbreviature: <YO777>

       _educ : integer;
        educ : array[1..2] of nstr;   // educts:  ME + YOU

       _prod : integer;               // quantite' des produits   (no. of products)
        prod : array[1..max_prod] of nstr; // products:  WE + THEY + US + ...
        stoi : array[1..max_prod] of real; // corr. stoichiometric coefficients:  1, many, some

        phys : string;                // right-part: {%DrSmSl} super*fast*rate;

        qatm_educ, qatm_prod : real;  // q-ty of atoms on both sides

        isrc : boolean;               // source-specification flag     (pipe me!)
        nsrc : integer;               // specification list entry no.

        iiex : boolean;               // is isotope exchange r-n?
        niex : integer;               // iex internal number

        itag : boolean;               // is tagged?
        etag : integer;               // is ever tagged? >0 tells how many times
        ntag : integer;               // tagging internal number

        if2t : boolean;               // are both educts tagged?
        if2s : boolean;               // are the educts same? (i.e. self-reacting)
        end;
   _eqs : word;                       // quantite' des reactions


// --- f-n returns back the number in eqs array for a certain reac abbr
function no_eqs(abbr : nstr) : integer;
var i : integer;
begin
no_eqs:=0;
if (trim(abbr)='') then exit;
for i:=1 to _eqs do
    if (abbr=eqs[i].abbr) then
       begin no_eqs:=i; exit; end;
end;

// --- f-n returns back the number in eqs array for a certain reac abbr
function no_educ(eqno : integer; ceduc : nstr) : integer;
var i : integer;
begin
no_educ:=0;
if (trim(ceduc)='') then exit;
with eqs[eqno] do
     for i:=1 to _educ do
         if (educ[i]=ceduc) then
            begin no_educ:=i; exit; end;
end;

// --- f-n returns back the true if a spec in the spc[] array
function in_spc(spec : nstr) : boolean;      // = (if (no_tsl(spec)>0))
var i : integer;
begin
in_spc:=false;
if (trim(spec)='') then exit;
for i:=1 to _spc do
    if (spec=spc[i].spec) then
       begin in_spc:=true; exit; end;
end;

// --- no of the species in the spc[] array
function no_spc(spec : nstr) : integer;
var i : integer;
begin
no_spc:=0;
if (trim(spec)='') then exit;
for i:=1 to _spc do
    if (spec=spc[i].spec) then
       begin no_spc:=i; exit; end;
end;

// --- checks whether the species acts in selected mechanism
function is_usedspec(cs : string) : boolean;
var s, r : word;
begin
is_usedspec:=false;
if (trim(cs)='') then exit;
for s:=1 to _eqs do
    with eqs[s] do
     begin
     if ((cs=educ[1]) or (cs=educ[2])) then
        begin
        is_usedspec:=true;
        exit;
        end;
     for r:=1 to _prod do
         if (cs=prod[r]) then
            begin
            is_usedspec:=true;
            exit;
            end;
     end;
end;

// --- no of spc[] entry if species has asrs<>'' (e.g. was added via ADD:SPC with {&spec})
function no_asrs(spec : nstr) : integer;
var i : integer;
begin
no_asrs:=0;
if (spec='') then exit;
for i:=1 to _spc do
    if (spec=spc[i].asrs) then
       begin no_asrs:=i; exit; end;
end;

// --- check whether certain reaction has a KIE
function is_kie(abbr : nstr) : boolean;
var i : word;
begin
is_kie:=false;
for i:=1 to _kie do
    if (kie[i].imec) then
       if (abbr=kie[i].abbr) then
          begin is_kie:=true; break; end;
end;

// --- check whether species acts in KIE reaction(s)
function in_kie(spec : nstr) : boolean;
var i, j : word;
begin
in_kie:=false;
for j:=1 to _kie do
    if (kie[j].imec) then
       with eqs[kie[j].eqno] do
            begin
            for i:=1 to 2 do
                if (in_tsl(educ[i]) and (educ[i]=spec)) then
                   begin in_kie:=true; exit; end;
            for i:=1 to _prod do
                if (in_tsl(prod[i]) and (prod[i]=spec)) then
                   begin in_kie:=true; exit; end;
            end;
end;

// --- la fonction, qui de'finit la condition de la re'action
// --- bool function determines criteria of the equation string recognition
function imcom_is_eqn(s : string) : boolean;
begin
imcom_is_eqn:=false;
if not(pos('//',trim(s))=1) then       // filtering commented lines
   if not(pos('!',trim(s))=1) then
      if (pos('<',s)<pos('>',s)) then  // some distinct markers
         if (pos('{%',s)>0) then
            imcom_is_eqn:=true;
end;

// --- determine whether criteria of the species string recognition
function imcom_is_spc(s : string) : boolean;
begin
imcom_is_spc:=false;
if not(pos('{',trim(s))=1) then        // filtering commented lines
   if not(pos('#',trim(s))=1) then     //   and inline keys
      if (pos('=',s)<pos(';',s)) then  //   some distinct markers
{        if (pos('{@',s)>0) then }     //   (not required)
            imcom_is_spc:=true;
end;



// --- reading and parsing equation file + checks for TSL given in CFG
procedure imcom_read_spc(fname : string);
var i, j, k, l : longint;
    a, s, siecntsl : ansistring;
    err : boolean;
    f : text;
begin

writeln('species file: ',fname);  // a ext

// error flag
err:=false;

imcom_check_files_exist([fname]);

// initializing
fillchar(spcfile,sizeof(eqnfile),0);
fillchar(spc,sizeof(spc),0);
_spc:=0; _spcfile:=0;

// species of interesting element in composition, but not in tsl
siecntsl:='';

// reading spc file by lines
assign(f,fname);
reset(f);
while not(eof(f)) do
      begin
      readln(f,s);

      inc(_spcfile);
      spcfile[_spcfile].line:=s;

      // if there were species in cfg to include, adding them to the top (foobar kpp)
      if (_addspc>0) then
         if (pos('#DEFVAR',upcase(trim(s)))=1) then
            begin
            inc(_spcfile);
            spcfile[_spcfile].line:='{------ [ '+tagname+' ] - additional species ------------------------------------}';

            for i:=1 to _addspc do
                begin
                inc(_spcfile);
                spcfile[_spcfile].line:=addspc[i];
                end;
            // make sure we do not add twice
            _addspc:=-1;
            end;
      end;
// closing eqn file
close(f);

// parsing line-by-line, searching for TSLs, plus pot. missed species
for l:=1 to _spcfile do
      begin
      show_runner(20);  // shows that we're not stuck parsing huge spc file

      s:=spcfile[l].line;

      if imcom_is_spc(s) then
         begin
         inc(_spc);
         spcfile[_spcfile].ispc:=true;
         spcfile[_spcfile].nspc:=_spc;

// - processing species:
         with spc[_spc] do

              begin

              // added spec. refspec. (if <>'', indicates trac. def. needs to be added)
              asrs:=trim(imcom_ext4marks(s,'{'+_asrssy,'}'));
              if (asrs<>'') then // cutting out (to prevent double detection in possible next cfg.)
                 begin
                 spcfile[l].line:=trimright(imcom_rem4marks(s,'{'+_asrssy,'}'));
                 s:=spcfile[l].line;
                 end;

              // getting caption for the regular
              capt[0]:=imcom_ext4marks(s,'{@','}');

              // couperons les commentaires dans les {}
              // cutting comments in curlies
              while (pos('{',s)>0) do
                    delete(s,pos('{',s),max(0,pos('}',s)-pos('{',s))+1);

              spec:=trim(copy(s,1,pos('=',s)-1));

              // getting composition for the regular
              s:=imcom_ext4marks(s,'=',';');
              comp[0]:=s;
              s:=delchars(s,' ')+'+';            // a marker for the last element detection

              // parsing composition
              while (s<>'') do
                    begin

                    // getting the qty number before the element
                    a:='';
                    while (s[1] in ['0'..'9']) do
                          begin
                          a+=s[1];
                          delete(s,1,1);
                          end;
                    j:=strtointdef(a,1);  // supposing no number = 1

                    // getting the element
                    a:='';
                    while not(s[1]='+') do
                          begin
                          a+=s[1];
                          delete(s,1,1);
                          end;
                    delete(s,1,1);

                    // calculating regular mol.mass
                    err:=true;
                    for i:=1 to _elems do
                        if (upcase(a)=upcase(elems[i].name)) then
                           begin
                           mass:=mass+j*elems[i].mass;
                           err:=false;
                           // checking if the species has the interesting element, but not in TSL
                           if not(in_tsl(spec)) and (i=_isoelem) then
                                 siecntsl+=spec+' ';
                           end;

                    if (err) then
                       writeln(' <!> info: for species '+spec+' element '+a+' is not known. this may affect molecule mass, => RMA');      // imcom_error

                    // linking SPC record to the TSL record
                    tslno:=no_tsl(spec);
                    if (tslno>0) then tsl[tslno].nspc:=_spc;

                    // composing isotopologues
                    for i:=1 to _isos do
                        begin
                        comp[i]+=' + ';

                        if (upcase(a)=upcase(isoelem)) then
                           begin
                           qatm:=j;                            // storing number of isotope atoms (for checking later)

                           if (j>1) then                          // multi-atom substitution
                              if (i>1) then
                                 begin
                                 comp[i]+=clsname[i]+isoelem+' + ';      // mixture of minor and major(s)
                                 if ((j-1)>1) then comp[i]+=inttostr(j-1);
                                 comp[i]+=clsname[1]+isoelem;
                                 end
                              else
                                  begin
                                  if (j>1) then comp[i]+=inttostr(j);  // only major(s)
                                  comp[i]+=clsname[1]+isoelem;
                                  end
                           else                                   // one-atom molecule substitution
                               comp[i]+=clsname[i]+isoelem;

                           // constructing caption: atom with superscripted class
                           //capt[i]:='^<'+clsname[i]+isoelem+'>';
                           end
                        else
                            begin                                 // not a substituted atom
                            if (j>1) then comp[i]+=inttostr(j);  // atom count
                            comp[i]+=a;
                            end;

                        // adding atom to the composition
{                        icapt[i]+=upcase(a);
                        if (j>1) then
                           icapt[i]+='_<'+inttostr(j)+'>';}

                        end;

                    end;

              for i:=1 to _isos do
                  begin
                  delete(comp[i],1,3);          // cut first " + "
                  capt[i]:=capt[0];
                  // trying to supstitute isotope
                  insert('^<'+clsname[i]+'>',capt[i],max(pos(isoelem,capt[0]),1));
                  // mecca spc-style caption brace
                  capt[i]:='{@'+capt[i]+'}';
                  end;

              end;
         end;

      end; // over spcfile

// un peu d'info
// a bit of information
writeln(' #',_spc,' species found');
writeln(' #',_tsl,' total in TSL');
writeln;
if (siecntsl<>'') then
   writeln(' <!> info: following species contain ',isoelem,' element, but not found in TSL: ',siecntsl,_LF);

// atom quantities checks
siecntsl:='';
for i:=1 to _tsl do
    if not(tsl[i].qatm>0) then  // assigning atoms count from SPC if not defined before through cfg
       begin
       tsl[i].qatm:=spc[tsl[i].nspc].qatm;
       siecntsl+=tsl[i].spec+'('+inttostr(tsl[i].qatm)+isoelem+') ';
       end;
if (siecntsl<>'') then
   writeln(' species that have atomic content picked from SPC: ',siecntsl,_LF);

siecntsl:='';
for i:=1 to _tsl do
    if not(tsl[i].qatm>0) then  // assigning atoms count = 1 for species not found both in SPC and TSL
       begin
       tsl[i].qatm:=1;
       spc[tsl[i].nspc].qatm:=1;
       siecntsl+=tsl[i].spec+' ';
       end;
if (siecntsl<>'') then
   writeln(' <!> warning: species with atomic content assigned to unity (not found in TSL/SPC): ',siecntsl,_LF);

siecntsl:='';
for i:=1 to _tsl do
    if (spc[tsl[i].nspc].qatm<>tsl[i].qatm) then
       siecntsl+=tsl[i].spec+' ('+inttostr(tsl[i].qatm)+'/'+inttostr(spc[tsl[i].nspc].qatm)+') ';
if (siecntsl<>'') then
   writeln(' <!> warning: following species have different ',isoelem,' content [using from TSL] (TSL/SPC) : ',siecntsl,_LF);

// masses check
siecntsl:='';
for i:=1 to _spc do
    if ((spc[i].mass < 1.0) or (spc[i].mass > 500.0)) then
       siecntsl+=spc[i].spec+' ('+floattostr(spc[i].mass)+') ';
if (siecntsl<>'') then
  writeln(' following species have strange mass [may be normal due to IGNORE]: ',siecntsl,' ');

// calculating masses of isotopologues (e.g. for RMA)
for i:=1 to _tsl do
    begin
    tsl[i].mass[0]:=spc[tsl[i].nspc].mass; // regular
    for k:=1 to _isos do                   // isotopologues
        tsl[i].mass[k]:=tsl[i].mass[0] - tsl[i].qatm*elems[_isoelem].mass  // removing qatm of weghted masses
                                       + (tsl[i].qatm-1)*isomass[1]        // adding (qatm-1) of lightest (abundant)
                                       + 1*isomass[k];                     // and one heavier (rare) mass
    end;

writeln;
writeln('read_spc: done');

writebreak;

end;

// returns the no. of product entry for species spc from eqs[i]
function no_prod(i : integer; spec : nstr) : integer;
var j : integer;
begin
with eqs[i] do
     for j:=1 to _prod do
         if (spec=prod[j]) then
            begin no_prod:=j; exit; end;
no_prod:=0;
end;

// --- lire et intepre'ter le fichier des equations
// --- reading and parsing equation file
procedure imcom_read_eqs(fname : string);

var i, j, k, l : longint;
    r, rm0, rmS, m1, m2 : real;
    s, a, ipf, ign_tag_eqs : ansistring;
    f : text;
    psign : integer;
    ex_prec, ex_prod : boolean;    // tagged-species exist on left/right side flags
    err, atag, fss, skp : boolean; // error flag, "already tagged" flag, source specification flag, on the skip-list

begin

//writeln('fichier des equations: ',fname);  // a ext
writeln('equations file: ',fname);  // a ext

// error flag
err:=false;

imcom_check_files_exist([fname]);

// lisons de l'info des reactions
// reading reactions info
assign(f,fname);
reset(f);

fillchar(eqnfile,sizeof(eqnfile),0);
fillchar(eqs,sizeof(eqs),0);
_eqs:=0; _eqnfile:=0; nooftagreac:=0;

// reading eqn file
while not(eof(f)) do
      begin
      readln(f,s);

      inc(_eqnfile);
      eqnfile[_eqnfile].line:=s;

      // if there were equations in cfg to include, adding them to the top (foobar kpp)
      if (pos('#EQUATIONS',upcase(trim(s)))=1) then
         if (_addeqs>0) then
            begin
            inc(_eqnfile);
            eqnfile[_eqnfile].line:='{------ [ '+tagname+' ] - additional reactions ----------------------------------}';

            for i:=1 to _addeqs do
                begin
                inc(_eqnfile);
                eqnfile[_eqnfile].line:=addeqs[i];
                end;
            // make sure we do not add twice
            _addeqs:=-1;
            end;
      end;
// closing eqn file
close(f);

// parsing line-by-line
for l:=1 to _eqnfile do
    begin
    show_runner(10);  // shows that we're not stuck reading huge eqn file

    s:=eqnfile[l].line;
    if imcom_is_eqn(s) then
       begin

       inc(_eqs);
       eqnfile[l].iseq:=true;
       eqnfile[l].eqno:=_eqs;
       eqs[_eqs].eqnf:=l;  // back-reference to eqnfile[] array

// - processing equation:
       with eqs[_eqs] do
            begin
            abbr:=imcom_ext4marks(s,'<','>');
            phys:=copy(s,pos(':',s)+1,length(s)-pos(s,':'));

            // "ever tagged" number: stored in the intermediate and final eqn by imtag
            etag:=strtointdef(imcom_ext4marks(s,'{'+_etagspsy,'}'),0);

            // couperons les commentaires dans les {}
            // cutting comments in curlies
            while (pos('{',s)>0) do
                  delete(s,pos('{',s),max(0,pos('}',s)-pos('{',s))+1);

            // identifying educt(s)
            if ((pos('+',s)<pos('=',s)) and (pos('+',s)>0)) then
               begin
               educ[1]:=trim(imcom_ext4marks(s,'>','+'));
               educ[2]:=trim(imcom_ext4marks(s,'+','='));
               _educ:=2;
               end
            else
                begin
                educ[1]:=trim(imcom_ext4marks(s,'>','='));
                educ[2]:='';
                _educ:=1;
                end;

            // couperons les co^te's gaushe
            // cutting left side
            delete(s,1,pos('=',s));
            // et droit
            // and right
            delete(s,pos(':',s),length(s)-pos(':',s)+1);

            s:=trim(s);
            if (length(s)>0) then
               if not(s[1] in ['+','-']) then s:='+'+s;  // adding a default '+' to the first product

            //writeln('}',s);

            i:=1;
            _prod:=0;
            fillchar(prod,sizeof(prod),0);
            while (i<length(s)) do
                  begin
                  inc(i);
                  case s[i-1] of      // last passed "+" or "-" belongs to the current product
                       '+' : psign:=1;
                       '-' : psign:=-1;
                  end; {case}

                  if (not(s[i] in ['+','-'])) then
                     prod[_prod+1]+=s[i];

                  if ((s[i] in ['+','-']) or (i=length(s))) then
                     begin
                     inc(_prod);
                     delete(s,1,length(prod[_prod])+1); s:=trim(s);
                     prod[_prod]:=trim(prod[_prod]);
                     //writeln('>',prod[_prod]);
                     //writeln(']',s);

                     // if the product has a stoi coeff
                     if (pos(' ',prod[_prod])>0) then
                        begin
                        stoi[_prod]:=strtofloat(copy(prod[_prod],1,pos(' ',prod[_prod])-1));
                        delete(prod[_prod],1,pos(' ',prod[_prod]));
                        end
                     else
                         stoi[_prod]:=1;

                     stoi[_prod]:=stoi[_prod]*psign;
                     i:=1;
                     end;

                   end;

            // la specification
            // source specification information
            eqs[_eqs].isrc:=false;
            for j:=1 to _src do
             // if (src[j].abbr=eqs[_eqs].abbr) then
                if pattmatch(eqs[_eqs].abbr,src[j].abbr,true) then  // matching reaction? (abbr. are case-sens.)
                   begin
                   // flag as "may have/has a source specification"
                   eqs[_eqs].isrc:=true;
                   eqs[_eqs].nsrc:=j;
                   src[j].eqno:=_eqs;
                   end;

            itag:=false;        // determined further
            ntag:=0;

            iiex:=false;      // not an isotope exchange r-n
            niex:=0;

            end;

// - EO initial processing

       end
    else
        eqnfile[l].iseq:=false;
    end;


// additional processing
// + un peu d'info:
ign_tag_eqs:='';

writeln('#',_eqs,' equations (showing only those containing species to '+_opmode+'):');
for i:=1 to _eqs do
    with eqs[i] do
         begin
         show_runner(10);  // shows that we're not stuck parsing huge eqn file

         // checking if the reaction in the eqn is already a taging r-n
         atag:=false;
         for j:=1 to _eqs do   // if among all of eqs. tags found one that matches current within trsyntax => tagged
             if (pos( substr(trsyntax,'@',eqs[j].abbr), eqs[i].abbr )>0 ) then
                begin
                atag:=true;
                ign_tag_eqs+=abbr+' ';
                break;
                end;
         if (atag) then continue;

         ex_prec:=false; ex_prod:=false;
         qatm_educ:=0; qatm_prod:=0;

         s:=format('%6s',[abbr])+': '; // write(abbr,': ');

         if (in_tsl(educ[1])) then
            begin
            s+=educ[1]+' '; //write(educ[1],' ');
            ex_prec:=true;
            qatm_educ+=tsl[no_tsl(educ[1])].qatm;
            end;
         if (in_tsl(educ[2])) then
            begin
            s+=educ[2]+' '; //write(educ[2],' ');
            ex_prec:=true;
            qatm_educ+=tsl[no_tsl(educ[2])].qatm;
            end;

         s:=s+'-> '; //write('-> ');

         for j:=1 to _prod do
             if (in_tsl(prod[j])) then
                begin
                if stoi[j]<>1 then s:=s+floattostr(stoi[j])+'*'; // write(floattostr(stoi[j]),'*');
                s+=prod[j]+' '; //write(prod[j],' ');
                ex_prod:=true;
                qatm_prod+=stoi[j]*tsl[no_tsl(prod[j])].qatm;
                end;

         // si l'equation n'a pas les molecules marques - passons
         if (not(ex_prec) and not(ex_prod)) then
            continue;

         // r-n has species from TSL, output info
         writeln('#',i:3,' @',etag:1,': ',trim(s));
         ipf:='       | <!> '; // additional info string prefix

         // checking whether the eqn in the skip-eqn-list
         skp:=false;
         for j:=1 to _skel do
             if pattmatch(abbr,skel[j],true) then  // matching reaction? (abbr. are case-sens.)
                begin
                writeln(ipf,'info: equation ',abbr,' has species to '+_opmode+', but is on the skip-list /',skel[j],'/ (may be OK)');
                skp:=true;
                continue;
                end;
         if (skp) then continue;

         // reaction is tagged:
         inc(nooftagreac);
         ntag:=nooftagreac;
         itag:=true;          // set "is tagged" flag
         etag:=etag+1;        // set "is ever tagged" flag, increase tot. number of tagging times
                              //  (relies on the hope that the PT abbreviation is equal for different configurations!)

         if2t:=(in_tsl(educ[1]) and in_tsl(educ[2]));  // both educts are tagged
         if2s:=(educ[1]=educ[2]);                      // self-reacting educts

         // checking budgeting
         for j:=1 to 2 do
             if (in_busl(educ[j])) then busl[no_busl(educ[j])].iloss:=true;
         for j:=1 to _prod do
             if (in_busl(prod[j])) then busl[no_busl(prod[j])].iprod:=true;


         // checking source specification
         if (isrc) then
            begin

            s:='';
            with src[eqs[i].nsrc] do
                 for j:=1 to _prod do                   // checking all products
                     if (in_tsl(prod[j])) then
                        begin
                        fss:=false;                     // flags whether a source is found for this product
                        for k:=1 to _trans do           // 160217: potential bugfix: used _trans-1 (with comment: -1 since we've created a new source)
                            if ( no_trans_dst(eqs[i].nsrc,k,prod[j])>0 ) and (trans[k].ib=0)  then   // there is destination, not an IB one
                               begin fss:=true; break; end;
                        if not(fss) then
                           if (src[eqs[i].nsrc].def>0) then        // if there is default source specified
                              with trans[src[eqs[i].nsrc].def] do
                                   begin
                                   inc(_dst);                      // adding transfer from the def. source
                                   dst[_dst].spec:=prod[j];
                                   dst[_dst].stoi:=1.0;            // so far, corrected further if necessary
                                   end
                           else
                               s+=prod[j]+' ';                     // keeping the list of products without source
                        end;

            // adding transfers for tagged educts without destination (important for destruction in tag without correction!)
            a:='';
            for j:=1 to _educ do
                if in_tsl(educ[j]) then
                   begin

                   k:=no_src_trans(eqs[i].nsrc,educ[j]);
                   if (k=0) then
                      with src[eqs[i].nsrc] do
                           begin
                           inc(_trans);
                           trans[_trans].src:=educ[j];
                           trans[_trans].weight:=1;
                           trans[_trans]._dst:=0;
                           a+=' '+educ[j];
                           end;
                   end;
            if (a<>'') then  writeln(ipf,'info: in ',abbr,': adding missing destruction term for:',a);

            if (s<>'') then  // there are unassigned products
               begin
               writeln(ipf,'error: in ',abbr,': no educt info for '+s+'(check eqn and configuration sources specification) -- reaction may be tagged improperly');
               //eqs[i].itag:=false; continue;
               //err:=true;
               end
            else
                // checking source stoichiometry (i.e. that partitioning factors of the source between destinations does not exceed unity)
                with src[eqs[i].nsrc] do
                     begin
                     s:='';
                     fss:=true;
                     while (fss) do
                           begin

                           // finding destination species not yet checked
                           fss:=false;
                           for j:=1 to _trans do
                             if (trans[j].ib=0) then      // skipping IB records
                               begin
                               for k:=1 to trans[j]._dst do
                                   if ( pos(trans[j].dst[k].spec,s)=0 ) then
                                      begin
                                      a:=trans[j].dst[k].spec;
                                      fss:=true;
                                      break;
                                      end;
                               if (fss) then break;
                               end;

                           if (fss) then
                              begin
                              r:=0; // qty. of molecules produced for the current destination (a)
                              for j:=1 to _trans do
                                  if (trans[j].ib=0) then    // skipping IB records
                                     if ( no_trans_dst(eqs[i].nsrc,j,a)>0 ) then
                                        r:=r+trans[j].dst[no_trans_dst(eqs[i].nsrc,j,a)].stoi;

                              // in [src] records coefficients are always molecular
                              // scaling excess over unity - as this "source" stoichiometry is later multiplied by
                              // (stochastic) branching and product molecular stoichiometry
                              // + values below unity are allowed to account for neglected destinations
                              if (r>1.0) then
                                 for j:=1 to _trans do
                                     if (trans[j].ib=0) then    // skipping IB records
                                        if ( no_trans_dst(eqs[i].nsrc,j,a)>0 ) then
                                           with trans[j].dst[no_trans_dst(eqs[i].nsrc,j,a)] do
                                                stoi:=stoi/r;

                              s+=' '+a; // remember processed spec

                              end;

                           end;

                     if (s<>'') then  writeln(ipf,'info: in ',abbr,': transfer stoichiometry is scaled to unity for:',s);
                     s:='';
                     end;

            // re-checking conservation of elemental/molecular transfer
            qatm_educ:=0; qatm_prod:=0;
            with src[eqs[i].nsrc] do
                 for j:=1 to _trans do
                     if (trans[j].ib=0) then    // skipping IB records
                        begin
                        qatm_educ+=tsl[no_tsl(trans[j].src)].qatm;
                        for k:=1 to trans[j]._dst do
                            begin
                            a:=trans[j].dst[k].spec;
                            qatm_prod+=trans[j].dst[k].stoi*( tsl[no_tsl(a)].qatm * eqs[i].stoi[no_prod(i,a)] );  // the latter term is from reg. r-n
                            end;
                        end;

            end; // if (isrc)


         // final checks

         // issue a warning if amount of consumed atoms is not equal to amount of produced
         if (abs(qatm_educ-qatm_prod)>eps_massbal) then
            begin
            write(ipf,'warning: in ',abbr,': ',isoelem,' element transfer is not conservative: ',
                  formatfloat('##0.###',qatm_educ), '->', formatfloat('##0.###',qatm_prod),' ');
            if (itransm=0) then write('(may be OK if molecular tagging is used)')
                           else write('(unlikely OK since isotope tagging is used)');
            writeln;
            writeln(ipf,'original r-n: ',eqnfile[eqs[i].eqnf].line);
            end;

         // issue a warning if a species comes from/to nothing
         if (ex_prec=not(ex_prod)) then
            if (ex_prec) then
               writeln(ipf,'warning: in ',abbr,': no product info (may be just a destruction)')
            else
                if not(eqs[i].isrc) then    // without the source specification? o-la-la
                   begin
{$IFNDEF IGNORE_NOSRC}
                   writeln(ipf,'error: in ',abbr,': no educt info (no source specification) -- cannot tag this reaction');
                   eqs[i].itag:=false; continue;
//                   err:=true;
{$ENDIF}
                   end
                else
                    writeln(ipf,'warning: in ',abbr,': no educt info from eqn (though specified in the configuration)');

         end;




j:=0; for i:=1 to _eqs do if (eqs[i].itag) then inc(j);
writeln('number of equations to tag: ',j);

{$IFDEF USE_DKRATE}
// is rate optimisatino is used, bubblesort comparison & correction of eqs. abbr (they need to be unique)
l:=0; s:='';
for i:=1 to _eqs-1 do
    begin

    // most time-taking part, showing runner
    show_runner(32);

    for j:=i+1 to _eqs do
        if (eqs[j].itag and eqs[i].itag) then
           if (upcase(eqs[j].abbr) = upcase(eqs[i].abbr)) then
              begin
              inc(l);
              eqs[j].abbr:=eqs[j].abbr+'_tag_'+inttostr(j);
              s+=eqs[j].abbr+' ';
              end;
    end;
write('corrected non-unique <abbr> reactions: ');
if (l=0) then writeln('none')
         else writeln('(',l,') '+s);
{$ENDIF}

writeln;

// isotope exchange reactions in the current mech
write('isotope exchange reactions (as for selected mechanism): ');
s:='';
// finding th last eq, -> l
l:=1; while (not(eqnfile[l].eqno=_eqs) and (l<=_eqnfile)) do inc(l); inc(l);
for i:=1 to _iex do
    with iex[i] do
         if (is_usedspec(tsl[exspec[1]].spec) and is_usedspec(tsl[exspec[2]].spec)) then
            begin
            // inserting a new line to the reactions list
            for j:=_eqnfile downto l do
                eqnfile[j+1]:=eqnfile[j];
            eqnfile[l].line:='';
            inc(_eqnfile);

            // adding new reaction
            inc(_eqs);
            eqs[_eqs].abbr:=abbr;
            eqs[_eqs].educ[1]:=tsl[exspec[1]].spec;
            eqs[_eqs].educ[2]:=tsl[exspec[2]].spec;
            eqs[_eqs]._prod:=2;
            eqs[_eqs].prod[1]:=tsl[exspec[2]].spec;
            eqs[_eqs].prod[2]:=tsl[exspec[1]].spec;
            eqs[_eqs].stoi[1]:=1.0;
            eqs[_eqs].stoi[2]:=1.0;

            eqs[_eqs].phys:=phys;

            eqs[_eqs].isrc:=false;
            eqs[_eqs].nsrc:=0;

            eqs[_eqs].iiex:=true;
            eqs[_eqs].niex:=i;

            inc(nooftagreac);
            eqs[_eqs].itag:=true;
            eqs[_eqs].etag:=1;    // since the IEX appears once (expected), we do not include it in the tag_k* stuff
            eqs[_eqs].ntag:=nooftagreac;

            // reactions file link
            eqnfile[l].line:='{ isotope exchange r-n '+abbr+' is to appear here }';
            eqnfile[l].iseq:=true;
            eqnfile[l].eqno:=_eqs;

            s+=abbr+' ';
            imec:=true;

            inc(l); // prepare for the next one
            end;
if (s='') then s:='none';
writeln(s);
writeln;

// making IESL
_iesl:=0; s:=''; a:='';
for i:=1 to _iex do
    if (iex[i].imec) then
       for j:=1 to 2 do
           if (pos('>'+tsl[iex[i].exspec[j]].spec+'<',s)=0) then
              begin
              inc(_iesl);
              iesl[_iesl]:=iex[i].exspec[j];
              s+='>'+tsl[iex[i].exspec[j]].spec+'<';
              a+=tsl[iex[i].exspec[j]].spec+' ';
              end;
write('species experiencing isotope exchange in selected mechanism: ');
if (_iesl>0) then writeln('(#',_iesl,') ',a) else writeln('none');

// source specification applied to the current reactions list
writeln('source specification (as for selected mechanism):');
for i:=1 to _eqs do
    if (eqs[i].isrc) and (eqs[i].itag) then
       with src[eqs[i].nsrc] do
            begin
            write(eqs[i].abbr:6,': ');

            for j:=1 to _trans do
              with trans[j] do
                begin
                if (ib>0) then write('(IB) ');       // indicate isotope branching
                if (weight<>1.0) then
                    write(formatfloat('##0.#####',weight)+' ');
                write(src,' ->');
                r:=0;
                for k:=1 to _dst do
                    begin
                    write(' ');
                    if (dst[k].stoi<>1.0) then
                       write(formatfloat('##0.#####',dst[k].stoi)+' ');
                    write(dst[k].spec);
                    r+=dst[k].stoi;
                    end;
                if (ib>0) and (abs(r-1.0)>eps_massbal) then a+=' '+src+'@'+eqs[i].abbr+'('+floattostr(r)+')';
                if (j<_trans) then write(', ');
                end;
            writeln;
            end;
if (a<>'') then writeln(' <!> warning: following sources have non-conservative IB:',a);
writeln;

// rejected "tagged" reactions list
if (length(ign_tag_eqs)>0) then
   writeln('ignored (already) tagged equations: ',ign_tag_eqs,_LF);

// making UTSL & fixed spec. count
_utsl:=0; s:='';
for i:=1 to _tsl do
    if (is_usedspec(tsl[i].spec)) then
       begin
       inc(_utsl);
       utsl[_utsl]:=i;
       s:=s+tsl[i].spec+' ';
       for j:=1 to _isos do
           if tsl[i].ifix[j] then 
              begin 
              inc(_fixs);
              break;
              end;
       end;
write('species to tag within the selected mechanism: ');
if (_utsl=0) then
  writeln('none   <!> warning: configuration '+cfgname+' has no species to tag!')
else
  writeln('(#',_utsl,') ',s);
writeln;

// indicating KIEs which are present in the current mechanism
s:='';
// regular reactions
for k:=1 to _kierec do
    for i:=1 to _eqs do
        if (eqs[i].itag) then                                        // checking only eqs. to tag
         if pattmatch(eqs[i].abbr,kierec[k].abbr,true) then           // matching reaction? (abbr. are case-sens.)
          if not(eqs[i].isrc) then                                     // checking source substitution differently (cf. below)
           begin
           for l:=1 to 2 do                                          // checking regular educts
               for j:=1 to _isos do
                   if pattmatch(tsl[no_tsl(eqs[i].educ[l])].isos[j],kierec[k].isot,false) then // regexp-matching isotopologue? (specs. are case-insens.)
                      begin
                      inc(_kie);
                      with kie[_kie] do
                           begin
                           abbr:=eqs[i].abbr;
                           isot:=tsl[no_tsl(eqs[i].educ[l])].isos[j];
                           expr:=kierec[k].expr;
                           eqno:=i;
                           imec:=true;
                           s:=s+abbr+'['+isot+'] ';
                           end;
                      end;
           end
        else { if (eqs[i].isrc) }
         with src[eqs[i].nsrc] do            // checking source substitution
              for l:=1 to _trans do                                  // e cycles through source specification
                  if (trans[l]._dst>0) then                          // avoiding sources without destination
                     for j:=1 to _isos do
//                       if (kie[k].isot=tsl[no_tsl(trans[l].src)].isos[j]) then
                         if pattmatch(tsl[no_tsl(trans[l].src)].isos[j],kierec[k].isot,false) then
                            begin
                            inc(_kie);
                            with kie[_kie] do
                                 begin
                                 abbr:=eqs[i].abbr;
                                 isot:=tsl[no_tsl(trans[l].src)].isos[j];
                                 expr:=kierec[k].expr;
                                 eqno:=i;
                                 imec:=true;
                                 s:=s+abbr+'['+isot+'] ';
                                 end;
                            end;
if (s='') then s:='none';
writeln('reactions with KIE within selected mechanism: ',s);

// checking reactions falling into the RMA selection
write('reactions with RMA KIE assigned: ');
for i:=1 to _eqs do
  if (eqs[i].itag) then                                               // checking only eqs. to tag
    for j:=1 to rmakie._include do
        if (pattmatch(eqs[i].abbr, rmakie.include[j], true)) then     // in case r-n abbrev. suits the wildcard given in RMA-include
           begin
           atag:=false;                                               // using this var as temporary boolean flag
           for k:=1 to rmakie._exclude do
               if (pattmatch(eqs[i].abbr, rmakie.exclude[k], true)) then // checking if r-n abbrev. in the exclude-list
                  atag:=true;
           for k:=1 to _kie do
               if (kie[k].abbr=eqs[i].abbr) then                      // anyway excluding reactions with already specified KIEs
                  atag:=true;
           if not(atag) then                                          // calculating and adding RMA KIE for the reaction
              begin
              write(eqs[i].abbr,' ');

              // adding RMA KIE to the regular kie list, he-he

              for l:=1 to 2 do
                  if (in_tsl(eqs[i].educ[l])) and                     // is an active spec?
                     not((eqs[i].educ[1]=eqs[i].educ[2]) and (l=2)) then   // avoiding doubling
                     for k:=2 to _isos do
                         begin

                         inc(_kie);
                         with kie[_kie] do
                              begin
                              abbr:=eqs[i].abbr;
                              eqno:=i;
                              imec:=true;
                              isot:=tsl[no_tsl(eqs[i].educ[l])].isos[k];

                              // picking up masses of isotopologues (or regular r-n partners)
                              m1:=tsl[no_tsl(eqs[i].educ[l])].mass[1];
                              if (in_tsl(eqs[i].educ[3-l])) then m2:=tsl[no_tsl(eqs[i].educ[3-l])].mass[1]  // using mol.mass with all lighter (abundant) isotopes
                                                            else m2:=spc[no_spc(eqs[i].educ[3-l])].mass;    // or with isotope-weighted mass for regular specs.

                              // change in mass from lighter to heavier (single substitution only!)
                              r:=-isomass[1]+isomass[k];

                              // checking division by 0
                              if ( ((m1  )*maxr(1,m2))<eps_zero ) or
                                 ( ((m1+r)*maxr(1,m2))<eps_zero ) then
                                 begin
                                 writeln('-> <!> warning: RMA KIE calculation for ',eqs[i].educ[l],'+',eqs[i].educ[3-l],
                                         ' (',tsl[no_tsl(eqs[i].educ[l])].isos[k],') most likely causes division by zero (calculating rm0/rmS); check species and cfg');
                                 rmS:=1; rm0:=1;
                                 end
                              else
                                  // calculating KIE using reduced masses (regular and substituted)
                                  begin
                                  rm0:=((m1  )+maxr(1,m2))/((m1  )*maxr(1,m2));
                                  rmS:=((m1+r)+maxr(1,m2))/((m1+r)*maxr(1,m2));
                                  if ( rm0 < eps_zero ) then
                                     begin
                                     writeln('-> <!> warning: RMA KIE calculation for ',eqs[i].educ[l],'+',eqs[i].educ[3-l],
                                         ' (',tsl[no_tsl(eqs[i].educ[l])].isos[k],') most likely causes division by zero (rm0 = 0); check species and cfg');
                                     rmS:=1; rm0:=1;
                                     end;
                                  end;
                              expr:='*'+floattostr(sqrt(rmS/rm0));
                              end;
                         end;

              end;
           end;
writeln;
writeln;


if (err) then
   imcom_error('read_eqs: errors detected, check the output and reactions file')
else
    writeln('read_eqs: done');

writebreak;

// some configuration is updated, updating replacements for templates
imcom_update_reps;

end;



// --- make_tracprop ----------------------------------------------------------
procedure imcom_make_tracprop(origTDname, addTDname : string;
                                        is_chemprop : boolean = false; // set to true when _chemprop (as opposed to _process) tbl input file is processed
                                          do_append : boolean = false; // set to true to append to the output file
                                           add_only : boolean = true); // set to fals to omit copying entries from input (only add new entries)
{$IFDEF TRACDEF_CHEMPROP}
// creates additional chemprop/process definition files used in MESSy TRACER
const tbl_delims : set of char = ['|'];
      tbl_remchr = '#'; // remark/comment char
      tbl_qtechr = '';  // quotes char
{$ELSE}
// creates additional tracers definition file to link with MESSy
const tbl_delims : set of char = ['&','\'];
      tbl_remchr = '%';  // remark/comment char
      tbl_qtechr = ''''; // quotes char
{$ENDIF}

const tbl_pad = ' '; // set to '' to remove space padding

// this f-n locates the no. of column with a given caption
function locate_column(src, caption : ansistring) : integer;
var w : integer;
begin
for w:=1 to wordcount(src, tbl_delims) do
  if ( trim(extractword(w,src,tbl_delims))=trim(caption) ) then
    begin locate_column:=w; exit; end;
locate_column:=0;
end;

// this function substitutes a certain column value with another
// style is similar to tracdef.tex, i.e.
// 'TRACER'  & col#2    & col#3  & etc.
function sub_column(src : ansistring; colno : integer; _repl : ansistring) : ansistring;
var p, n, w : integer;
    s : ansistring;
begin
(*
s:=src;
// searching for a column and getting its width
p:=0; w:=0; n:={$IFDEF TRACDEF_CHEMPROP}-1{$ELSE}0{$ENDIF};  // tables in chemprop start with a delimiter
while ((p<length(src)) and (n<colno)) do
      begin
      inc(p);                      // points to the end of col.
      inc(w);                      // curr. col width
      if (s[p] in tbl_delims) then
         begin
         inc(n);                   // passed column no.
         if (n<colno) then w:=0;   // resetting column width
         end;
      end;
if (n=colno) then                  // we found the col.!
   begin
   dec(w); dec(p,w);
   delete(s,p,w);                  // removing old cont., putting new
   insert(format('%-'+inttostr(w)+'s',[_repl]),s,p);
   end;
sub_column:=s;
*)
if (colno>0) then
  begin
  w:=length(extractword(colno,src,tbl_delims));
  p:=wordposition(colno,src,tbl_delims);
  delete(src,p,w);                                   // removing old
  insert(format('%-'+inttostr(w)+'s',[tbl_pad+_repl+tbl_pad]),src,p);  // insterting new
  end;
sub_column:=src;
end;

var f : text;
    tdfile : array of ansistring;  // careful here with long strings in *.tex !!!
   _tdfile : integer;              // # of lines in the file
    i, j, k, fn : integer;
    otrac, refspec : nstr;         // original tracer, and, if present, its refspec property (to drill down to refspec that is present in chemprop, only 1 level!)
    s, ssave : ansistring;
    rmm_fmt, rmm_fmt_full : string; // R_molarmass value format (rounded and full)

// now these are the indices in cols_adj_caps
const col_bname = 1; //= 1;         // basename
      col_cprop = 2;                // <CHEMPROP>
      col_med   = 3; //= 3;         // media
      col_vini  = 4; //= 4;         // vini
      col_rmm   = 5; //= 17;        // R_molarmass
      col_rspec = 6;                // refspec for TRACER
      col_sname = 7;                // subname

var cols_qtotal: integer;   // total # of the columns in .tex to process
    col_name : integer;     // actual entry for name (either _bname or _cprop)
    last_entry : integer;   // last detected tracer entry

 // columns detected in tracdef
var cols : array of string;
   _cols : integer;

 // columns on which we do adjustments (rest is copied or left out)
const cols_ADJ_caps : array[1..7] of string = ('basename','<CHEMPROP>','medium','R_vini','R_molarmass','refspec','subname'); //  'R_Henry_T0','R_Henry_Tdep','R_Sechenov','R_alpha_T0','R_alpha_Tdep','R_K_acid','R_K_acid2','R_Vmol_bp','R_psat','R_psat_Tdep','R_pss','R_dryreac_sf'
  var cols_ADJ_nos  : array[1..7] of integer;
{$IFDEF TRACDEF_CHEMPROP}
  var cols_ADJ_nos_def : array[1..7] of integer = (-1    , -1         , 5      , 6      , 3           , 4       , 2       );
{$ELSE}
  var cols_ADJ_nos_def : array[1..7] of integer = (-1    , -1         , 3      , 4      , 17          , -1      , -1      );
{$ENDIF}

 // columns for which PTs have to have OFF values
const cols_OFF_caps : array[1..10] of string = ('I_force_init','I_advect','I_convect','I_vdiff','I_wetdep','I_drydep','I_sedi','I_scav','I_mix','I_force_col'); //,'I_integrate');
  var cols_OFF_nos  : array[1..10] of integer;
{$IFDEF TRACDEF_CHEMPROP}
  var cols_OFF_nos_def : array[1..10] of integer = (     6    , 7        , 8         , 9       , 10       , 11       , 12     , 13     , 14    , 15  );
{$ELSE}
  var cols_OFF_nos_def : array[1..10] of integer = (     5    , 6        , 7         , 8       , 9        , 10       , 11     , 12     , 13    , 14  );
{$ENDIF}

var isout : array[1..max_tsl] of boolean; // is tracdef for a given species written out?
    H2O_present : boolean;

begin

// defines the format of R_molarmass values
rmm_fmt:='%-.'+inttostr(_roundmass)+'f'; // +'_dp' // yields rounded to _roundmass's digit 0.##0 format
rmm_fmt_full:='%-.8f';                   // yields 0.00000001 format

writeln('imcom_make_tracprop(',origTDname,',',addTDname,'):');
// reading original file to get samples from
imcom_check_files_exist([origTDname]);
assign(f,origTDname);

// reading in the file
reset(f);
_tdfile:=0; setlength(tdfile,0);
while not(eof(f)) do
  begin
  inc(_tdfile);
  setlength(tdfile,1+_tdfile);
  readln(f,tdfile[_tdfile]);
  end;
close(f);

// opening output file handle
assign(f,addTDname); if (do_append) then append(f) else rewrite(f);
// info
writeln(f,tbl_remchr,' this file was created/modified by [',{$IFDEF EMB}'embudget'{$ELSE}'imtag'{$ENDIF},']');
writeln(f,tbl_remchr,' configuration: ',cfgname);
writeln(f,tbl_remchr,' created: ',datetimetostr(now));
writeln(f,tbl_remchr,'');

// setting data columns to process and identifying their nos.
// all [ 'basename','subname','medium','vini','I_force_init','I_advect','I_convect','I_vdiff','I_wetdep','I_drydep','I_sedi','I_scav','I_mix','I_force_col','I_integrate','I_timefilter','I_method','S_model','I_mode',' I_sol' ]
cols_qtotal:=0; // total # of the columns in tex/tbl
for j:=1 to length(cols_ADJ_caps) do cols_ADJ_nos[j]:=0;
for j:=1 to length(cols_OFF_caps) do cols_OFF_nos[j]:=0;
for i:=1 to _tdfile do
  begin
  s:=tdfile[i];

  // columns with values to be adjusted
  for j:=1 to length(cols_adj_nos) do
    if (cols_ADJ_nos[j]=0) then cols_adj_nos[j]:=locate_column(s,cols_ADJ_caps[j]);

  // columns with values to be set OFF
  for j:=1 to length(cols_OFF_nos) do
    if (cols_OFF_nos[j]=0) then cols_OFF_nos[j]:=locate_column(s,cols_OFF_caps[j]);

  // max. columns no. met in file
  cols_qtotal:=max(cols_qtotal,wordcount(s,tbl_delims));
  end;

// some checks
if not(cols_qtotal>1) then
  imcom_error('imcom_make_tracprop: could not detect any column')
else
  begin
  s:='adjust(';
  for j:=1 to length(cols_ADJ_nos) do
    if (cols_adj_nos[j]<>0) then s+=cols_ADJ_caps[j]+'['+inttostr(cols_ADJ_nos[j])+'] ';
  s:=trim(s)+'), off(';
  for j:=1 to length(cols_OFF_nos) do
    if (cols_OFF_nos[j]<>0) then s+=cols_OFF_caps[j]+'['+inttostr(cols_OFF_nos[j])+'] ';
  s:=trim(s)+')';
  writeln('- detected columns: ',s);
end;

// required
col_name:=0;
if (cols_adj_nos[col_cprop]>0) then col_name:=col_cprop;
if (cols_adj_nos[col_bname]>0) then col_name:=col_bname;
if (col_name=0) then
  imcom_error('imcom_make_tracprop: could not detect either of <CHEMPROP> or basename column');
// if this is a process file, check for subname column
if not(is_chemprop) and (cols_adj_nos[col_sname]=0) then
  writeln(' <!> warning: could not detect subname column (may cause errors with `process` input)');

{$IFNDEF TRACDEF_CHEMPROP}
// attempt to substitute not-found column nos. with defaults
s:='';
for j:=1 to length(cols_ADJ_nos) do
  if (cols_ADJ_nos[j]=0) and (cols_ADJ_nos_def[j]>0) then
     begin cols_adj_nos[j]:=cols_adj_nos_def[j]; s+=' '+cols_ADJ_caps[j]+'['+inttostr(cols_ADJ_nos_def[j])+']'; end;
for j:=1 to length(cols_OFF_nos) do
  if (cols_OFF_nos[j]=0) and (cols_OFF_nos_def[j]>0) then
     begin cols_OFF_nos[j]:=cols_OFF_nos_def[j]; s+=' '+cols_OFF_caps[j]+'['+inttostr(cols_OFF_nos_def[j])+']'; end;
if (s<>'') then
   writeln(' <!> warning: could not detect all ADJ/OFF columns; assumed default positions for: ',s);
{$ENDIF}

// re-detecting all column names from the line containing <CHEMPROP> or basename
_cols:=0;
//writeln('searching for: ',cols_adj_caps[col_name],' cols_adj_nos[col_cprop]=',cols_adj_nos[col_cprop],' cols_adj_nos[col_bname]=',cols_adj_nos[col_bname]);
for i:=1 to _tdfile do
  if (locate_column(tdfile[i],cols_ADJ_caps[col_name])>0) then
    begin
    for j:=1 to wordcount(tdfile[i],tbl_delims) do
      begin
      inc(_cols); setlength(cols,1+_cols);
      cols[j]:=trim(extractword(j,tdfile[i],tbl_delims));
      end;
    break;
    end;
// checks
if (_cols<1) then
  imcom_error('imcom_make_tracprop: could not detect column captions (possible double entry of basename or <CHEMPROP>)');

(*
k:=0; // no of detected columns
for j:=3 to length(cols_adj_nos) do       // columns with values to be adjusted
  if (cols_adj_nos[j]>0) then inc(k);
if (k=3) then writeln('<!> imcom_make_tracprop/addtracdef: could not detect all columns to be adjusted');
k:=0;
for j:=1 to length(cols_OFF_nos) do       // columns with values to be set OFF
  if (cols_OFF_nos[j]>0) then inc(k);
if (k=0) then writeln('<!> imcom_make_tracprop/addtracdef: could not detect all columns to be set OFF');
*)

// detecting last table entry and H2O
last_entry:=0;
H2O_present:=false;
for i:=1 to _tdfile do
  begin
  s:=trim(imcom_rem4marks(extractword(cols_adj_nos[col_name],tdfile[i],tbl_delims),'{','}'));
  if (tbl_qtechr<>'') then s:=replacestr(s,tbl_qtechr,''); // remove quotes, if any
  if in_spc(s) then last_entry:=i;
  if (upcase(s)='H2O') then H2O_present:=true;
  end;
if (last_entry=0) then
   imcom_error('imcom_make_tracprop: could not detect last entry');
// done with column/entries identification

{$IFDEF TRACDEF_ADD_H2O}
// adding a 'fake' H2O tracer in order to create tagged tracers
if (not(H2O_present) and (in_tsl('H2O') or (no_asrs('H2O')>0))) then
  begin
  writeln(' <!> warning: imcom_make_tracprop: did not detect H2O (may be OK), adding');
  // inserting 2 new lines after the last_entry
  inc(_tdfile,2); setlength(tdfile,1+_tdfile);
  for i:=last_entry+3 to _tdfile do
    tdfile[i]:=tdfile[i-2];

  // info
  tdfile[last_entry+1]:=tbl_remchr+'--- this H2O tracer was added by imtag -';
  while (length(s)<length(tdfile[last_entry])) do s+='-';

  // using last_entry to compose an H2O
  s:=tdfile[last_entry];

  // removing all data (except the medium)
  for k:=1 to cols_qtotal do
    if (k<>cols_adj_nos[col_med]) then
      s:=sub_column(s, k, '');

  // setting the name
  s:=sub_column(s, cols_adj_nos[col_name], tbl_qtechr+'H2O'+tbl_qtechr);

  // changing R_molarmass (if isotope masses were specified)
  s:=sub_column(s, cols_adj_nos[col_rmm], {$IFDEF TRACDEF_CHEMPROP}'H2O'{$ELSE}'MH*2._dp+MO'{$ENDIF});

  // adjusting the I_ parameters
  //cols_OFF_caps : array[1..10] of string = ('I_force_init','I_advect','I_convect','I_vdiff','I_wetdep','I_drydep','I_sedi','I_scav','I_mix','I_force_col'); //,'I_integrate');
  for k:=2 to 4 do
    s:=sub_column(s, cols_OFF_nos[k], 'ON');
  for k:=5 to 8 do
    s:=sub_column(s, cols_OFF_nos[k], 'OFF');

  // adding
  inc(last_entry,2); tdfile[last_entry]:=s;
  end;
{$ENDIF}

// reset info on output tracdefs for specs
fillchar(isout,sizeof(isout),0);

ssave:='';

// processing sequence: finding a spec, if it is used in the tagging
// addtracdef: putting new tracers according to the classes into new file
//   tracprop: adding new tracers to the original file
for i:=1 to _tdfile do
    begin

    if (add_only) then         // retain header lines
      begin
      {$IFDEF TRACDEF_CHEMPROP}
      // output original line in tbl if this is a header
      for j:=1 to length(cols_adj_caps) do
        if (locate_column(tdfile[i],cols_adj_caps[j])>0) then
           begin writeln(f,tbl_remchr+rightstr(tdfile[i],length(tdfile[i])-1)); break; end;
      {$ELSE}
      // copying first 4 lines (supp. to be the table header)
      if (i<=min(4,_tdfile)) then writeln(f,tdfile[i]);
      {$ENDIF}
      end
    else
      writeln(f,tdfile[i]);    // retain original line

    // skipping commented lines
    if (copy(trim(tdfile[i]),1,1)=tbl_remchr) then continue;

    // drill down regular tracer name
    otrac:=trim(imcom_rem4marks(extractword(cols_adj_nos[col_name],tdfile[i],tbl_delims),'{','}'));
    if (tbl_qtechr<>'') then otrac:=replacestr(otrac,tbl_qtechr,''); // remove quotes, if any
    otrac:=trim(otrac);
    // and its refspec, if present
    refspec:=trim(imcom_rem4marks(extractword(cols_adj_nos[col_rspec],tdfile[i],tbl_delims),'{','}'));
    if (tbl_qtechr<>'') then refspec:=replacestr(refspec,tbl_qtechr,''); // remove quotes, if any
    // setting refspec to otrac, if it is nil
    if (refspec='') then refspec:=otrac;
    // if subname is present, add it
    if (cols_adj_nos[col_sname]<>0) then
      begin
      s:=trim(imcom_rem4marks(extractword(cols_adj_nos[col_sname],tdfile[i],tbl_delims),'{','}'));
      if (tbl_qtechr<>'') then s:=trim(replacestr(s,tbl_qtechr,'')); // remove quotes, if any
      if (s<>'') then otrac+='_'+s;
      end;

    // if specied is used in tagging/budgeting or added with reference, add tracer(s)
    if {$IFDEF EMB}(in_tsl(otrac)){$ELSE}(in_utsl(otrac)){$ENDIF} or (no_asrs(otrac)>0) then
      begin

      // saving the sample of a tracer def. line
      ssave:=tdfile[i];

      // line referring to the regular spec
      s:=tbl_remchr+'- '+otrac+' -';
      while (length(s)<length(tdfile[i])) do s+='-';
      writeln(f,s);
      // + commented original entry
      s:=ssave;
      if (s[1]<>tbl_qtechr) then s[1]:=tbl_remchr else insert(tbl_remchr,s,1);
      writeln(f,s);

      // custom added species (with [ADD:SPC], spc[].asrs<>''), refspec=otrac in case no refspec given
      if (no_asrs(refspec)>0) then
         begin
         for j:=1 to _spc do
             if (spc[j].asrs=otrac) then
                begin
                s:=ssave;
                // changing tracer name
                s:=sub_column(s, cols_adj_nos[col_name], ' '+tbl_qtechr+(spc[j].spec)+tbl_qtechr);
{$IFDEF TRACDEF_CHEMPROP}
                // removing all chem. prop. data I_, R_, S_
                for k:=1 to _cols do
                    if (is_chemprop and pattmatch(cols[k],'[IRS]_.+',true)) then // removing all properties (tracer will be created using refspec)
                       s:=sub_column(s, k, '');                                  // for _process tbl, retain everything
                // adding reference to the ref. species
                s:=sub_column(s, cols_adj_nos[col_rspec], tbl_qtechr+otrac+tbl_qtechr);
{$ELSE}
                // retain all data in case of tracdef (refspec not implemented)
//              for k:=1 to cols_qtotal do
//                  if (k<>cols_adj_nos[col_med]) and (k<>cols_adj_nos[col_rmm]) then
//                     s:=sub_column(s, k, ' ');
{$ENDIF}
                // write out
                writeln(f,s);
                end;
         // cycle if the added species is not in [u]tsl
         if not({$IFDEF EMB}(in_tsl(otrac)){$ELSE}(in_utsl(otrac)){$ENDIF}) then continue;
         end;

{$IFDEF EMB}

      // --- part compiled for embudget --------------------
      k:=no_tsl(otrac);
      isout[k]:=true; // flag as written out

      with tsl[k],tsla[k] do
        for j:=1 to _cats do
          begin

          s:=tdfile[i];

{$IFDEF TRACDEF_CHEMPROP}
          // removing all chem. prop. data I_, R_, S_
          for k:=1 to _cols do
            if (pattmatch(cols[k],'[IRS]_.+',true)) then
              s:=sub_column(s, k, '');
          // adding reference to the budgeted regular tracer
          s:=sub_column(s, cols_adj_nos[col_rspec], tbl_qtechr+refspec+tbl_qtechr);
{$ELSE}
          // removing all initial data (except the medium and R_molarmass)
          for k:=1 to cols_qtotal do
            if (k<>cols_adj_nos[col_med]) and (k<>cols_adj_nos[col_rmm]) then
              s:=sub_column(s, k, ' ');
{$ENDIF}

          // turning off all required I_ parameters
          for k:=low(cols_OFF_nos) to high(cols_OFF_nos) do
              s:=sub_column(s, cols_OFF_nos[k], 'OFF');

          if (iloss) and (uloss[j]) then
            begin
            // tracer name
            s:=sub_column(s, cols_adj_nos[col_name], ' '+tbl_qtechr+substr(ptlsyntax,'@',ncat[j]+spec)+tbl_qtechr);
            // write down
            writeln(f,s);
            end;
          if (iprod) and (uprod[j]) then
            begin
            // tracer name
            s:=sub_column(s, cols_adj_nos[col_name], ' '+tbl_qtechr+substr(ptpsyntax,'@',ncat[j]+spec)+tbl_qtechr);
            // write down
            writeln(f,s);
            end;
          end;

{$ELSE} // EMB

      // --- part compiled for imtag -----------------------
      k:=no_utsl(otrac);
      isout[k]:=true; // flag as written out

      // tagging counterparts / isotopologues
      for j:=1 to _isos do
        begin

        s:=tdfile[i];

{$IFDEF TRACDEF_CHEMPROP}
        // removing all chem. prop. data I_, R_, S_
        for k:=1 to _cols do
          if (pattmatch(cols[k],'[IRS]_.+',true)) then
            s:=sub_column(s, k, '');
        // adding reference to regular tracer
        s:=sub_column(s, cols_adj_nos[col_rspec], tbl_qtechr+refspec+tbl_qtechr);
{$ENDIF}

        // tracer name
        s:=sub_column(s, cols_adj_nos[col_name], ' '+tbl_qtechr+tsl[no_tsl(otrac)].isos[j]+tbl_qtechr);

        // removing initial conc. value
        s:=sub_column(s, cols_adj_nos[col_vini], '');

        // changing R_molarmass (if isotope masses were specified)
        if (isomass[j]>0.0) then
          s:=sub_column(s, cols_adj_nos[col_rmm], format( rmm_fmt, [ tsl[no_tsl(otrac)].mass[j] ] ) );

        writeln(f,s);
        end;
      end;  // in_[u]tsl


    // budgeting tracers (PTs)
    if (in_busl(otrac)) then
      with busl[no_busl(otrac)] do
        for j:=1 to _isos do
          begin

          s:=tdfile[i];

{$IFDEF TRACDEF_CHEMPROP}
          // removing all chem. prop. data I_, R_, S_
          for k:=1 to _cols do
            if (pattmatch(cols[k],'[IRS]_.+',true)) then
              s:=sub_column(s, k, '');
          // adding reference to the regular species
          s:=sub_column(s, cols_adj_nos[col_rspec], tbl_qtechr+refspec+tbl_qtechr);
{$ELSE}
          // removing all initial data (except the medium and R_molarmass)
          for k:=1 to cols_qtotal do
            if (k<>cols_adj_nos[col_med]) and (k<>cols_adj_nos[col_rmm]) then
              s:=sub_column(s, k, ' ');
{$ENDIF}

          // removing all chem. prop. data I_, R_, S_
          for k:=1 to _cols do
            if (pattmatch(cols[k],'[IRS]_.+',true)) then
              s:=sub_column(s, k, '');

          // removing initial conc. value
          s:=sub_column(s, cols_adj_nos[col_vini], '');

          // turning off all required I_ parameters
          for k:=low(cols_OFF_nos) to high(cols_OFF_nos) do
            s:=sub_column(s, cols_OFF_nos[k], 'OFF');

          if (iloss) then
            begin
            // loss PT tracer name
            s:=sub_column(s, cols_adj_nos[col_name], ' '+tbl_qtechr+substr(ptlsyntax,'@',clsname[j]+busl[no_busl(otrac)].spec)+tbl_qtechr);
            writeln(f,s);
          end;
          if (iprod) then
            begin
            // production PT tracer name
            s:=sub_column(s, cols_adj_nos[col_name], ' '+tbl_qtechr+substr(ptpsyntax,'@',clsname[j]+busl[no_busl(otrac)].spec)+tbl_qtechr);
            writeln(f,s);
            end;
          end;

    // other tracers (as last entries in table)
    if (i=last_entry) then
      begin

{$IFDEF USE_PT_UPL}
      // ----------- adding unaccounted production/loss tracers --------------
      s:=tbl_remchr+'--- '+cfgname+' unaccounted production/loss PTs ---';
      while (length(s)<length(ssave)) do s+='-';
      writeln(f,s);

      if (ssave='') then ssave:=tdfile[i];

      // unaccounted {%ATOM} tracers: turning off all processes and assigning mass of the resp. isotope (element)
      if (ssave<>'') then
        begin

        s:=ssave;

        // removing all initial data (except the medium)
        for k:=1 to cols_qtotal do
          if (k<>cols_adj_nos[col_med]) then
            s:=sub_column(s, k, '');

        // changing Rmolarmass    // maxr is to assure that there is no zero-mass
        s:=sub_column(s, cols_adj_nos[col_rmm], format( rmm_fmt_full, [ maxr(elems[-1].mass,elems[_isoelem].mass) ] ));

        // turning off all I_ parameters
        for k:=low(cols_OFF_nos) to high(cols_OFF_nos) do
          s:=sub_column(s, cols_OFF_nos[k], 'OFF');

        // PTLU tracer
        s:=sub_column(s, cols_adj_nos[col_name], tbl_qtechr+substr(ptlsyntax,'@','U'+cfgname+isoelem)+tbl_qtechr);
        writeln(f,s);

        // PTPU tracer
        s:=sub_column(s, cols_adj_nos[col_name], tbl_qtechr+substr(ptpsyntax,'@','U'+cfgname+isoelem)+tbl_qtechr);
        writeln(f,s);

        // and "isotopologues"
        for j:=1 to _isos do
          begin
          // changing Rmolarmass    // maxr is to assure that there is no zero-mass
          s:=sub_column(s, cols_adj_nos[col_rmm], format( rmm_fmt_full, [ maxr(elems[-1].mass,isomass[j]) ] ));
          // PTLU tracer
          s:=sub_column(s, cols_adj_nos[col_name], tbl_qtechr+substr(ptlsyntax,'@','U'+cfgname+clsname[j])+tbl_qtechr);
          writeln(f,s);
          // PTPU tracer
          s:=sub_column(s, cols_adj_nos[col_name], tbl_qtechr+substr(ptpsyntax,'@','U'+cfgname+clsname[j])+tbl_qtechr);
          writeln(f,s);
          end;
        end
      else
          writeln(' <!> warning: could not find any suitable tracer definition line in <'+
                    origTDname+'>; unaccounted prod/loss tracers are not listed in tbl/tracdef.');
{$ENDIF}

{$ENDIF} // EMB
      end;

    end; // for i:=1 to _tdfile do

// last row line (cosmetics)
s:=tbl_remchr+'--';
while (length(s)<length(ssave)) do s+='-';
writeln(f,s);

close(f);
// cleanup
_tdfile:=0; setlength(tdfile,_tdfile);
_cols:=0; setlength(cols,_cols);

// check whether tracdefs for all species were written out
s:='';
for k:=1 to {$IFDEF EMB}_tsl{$ELSE}_utsl{$ENDIF} do
  if not(isout[k]) then s+=' '+tsl[{$IFDEF EMB}k{$ELSE}utsl[k]{$ENDIF}].spec;
if (s<>'') then
  writeln(' <!> warning: imcom_make_tracprop/addtracdef: could not add entries for',s,' in <',origTDname,'>');

writeln('- done');

end;


// - util ---------------------------------------------------------------------

// -- function to wrap the lines, (un)conditionally by &
                           // auto &   conditional indents
function wrap(a : ansistring; aindent, cindent : word) : ansistring;
var s : ansistring;
    c, p : integer;
    w : boolean;
    lead, blank : string;
begin

s:='';
c:=1;
w:=false;

setlength(blank,aindent+cindent);
fillchar(blank[1],aindent+cindent,#32);
//a:=copy(blank,1,aindent)+trimleft(a);        // first line indent

while (c<length(a)) do
      begin

      p:=0;

      lead:=copy(blank,1,aindent);

      while ((c<length(a)) and not(w)) do
            begin
            if (a[c] in ['+','-',',']) then p:=c;       // linebreak allowed at "+" and ',' only
            if ((c>=wrapmargin) and (p>0)) then w:=true;
            // unconditional wrapping at ampersand
            if (a[c]='&') then
               begin
               delete(a,c,1);
               insert(' ',a,c);
               p:=c;
               w:=true;
               lead:=copy(blank,1,cindent);
               end;
            inc(c);
            end;

      if not(c<length(a)) then
         begin
         s:=s+a;
         break;
         end;

      c:=1; w:=false;

      s:=s+copy(a,1,p)+' &'+_LF;
      delete(a,1,p);
      a:=trimleft(a);
      while (a[1]='&') do
            delete(a,1,1);     // avoid exorbitant ampersands from code generator
      insert(lead,a,1);
      end;

wrap:=s;

end;


// - merged subroutines --------------------------------------------------------


// imcom_condition -------------------------------------------------------------
// returns a certain condition comparison result
function imcom_condition(c, v : string) : boolean;
var z, n : boolean;
    i : integer;
begin

z:=false;

if c[1]='!' then                                    // negated condition
   begin n:=true; delete(c,1,1); end
else n:=false;

c:=uppercase(c);
v:=trim(v);

if c='CASE' then z:=(v='TAG') else
if c='ATOM' then z:=pattmatch(isoelem,v,true) else  // element is
if c='CONF' then z:=pattmatch(cfgname,v,true) else  // configuration is ...
if c='TAG'  then z:=pattmatch(tagname,v,true) else  // tagging name
if c='ISO'  then z:=(_isos=strtointdef(v,-1)) else  // number of classes is equal to
if c='KIE'  then z:=(_kie>0) else                   // KIE present in the system
if c='SRC'  then z:=(_src>0) else                   // specification of sources is present
if c='REAC' then
   begin
   for i:=1 to _eqs do
       if eqs[i].itag then
          if pattmatch(eqs[i].abbr,v,true) then
             z:=true;
   end else
if c='SPEC' then z:=is_usedspec(v) else             // if in the used species list
if c='CPAR' then  // if an expression is in cfg parameters [CPARAMS] list
   begin
   z:=false;
   for i:=1 to wordcount(cparams,_delims) do
       if pattmatch(extractword(i,cparams,_delims),v,true) then
          z:=true;
   end else
   begin
   writeln(' <!> warning: imcom_condition( param: ',c,', value: ',v,' ) is not recognized');
   z:=false;
   end;

if (n) then z:=not(z);
imcom_condition:=z;

end;



// imcom_tagspecs (tagged species list) ----------------------------------------
// outputs the tagged species names list according to specimen with @, #/$ (class #/name)
function imcom_make_tagspecs_list(pref : string) : ansistring;
var i, k : word;
    a : ansistring;
begin

a:='       ';
// tracers names (@), no cycling through classes
if ((pos('#',pref)+pos('$',pref))=0) then
   for i:=1 to _utsl do
       a:=a+''''+format('%-'+inttostr(sizeof(nstr)-1)+'s',[substr(pref,'@',tsl[utsl[i]].spec)])+''', '
else // or cycling classes, if marker for no.(#)/class name($) is found
    for k:=1 to _isos do
        begin
        for i:=1 to _utsl do
            a:=a+''''+format('%-'+inttostr(sizeof(nstr)-1)+'s',
               [ substr(substr(substr(pref,'@',tsl[utsl[i]].spec),'#',inttostr(k)),'$',clsname[k]) ])+''', ';
        end;
setlength(a,length(a)-2);
a:=a+' &';

imcom_make_tagspecs_list:=wrap(a,7,7);
end;


// -----------------------------------------------------------------------------
// outputs an array of indices of species acting in reactions with KIE
// prefix: ind_d for doubling, ind_t for tagging
function imcom_make_kie_iex_relspecs(reftrac, tagtrac : string) : ansistring;
var i, j, k, r, s : word;
    a, b, c, out : string;
begin

out:='';

// KIE part

// making list of reactions with KIE in the mech
r:=0; s:=0; k:=0; a:=''; b:=''; c:='';
for i:=1 to _kie do
    if (kie[i].imec) then
       if (pos(kie[i].abbr,b)=0) then
          begin
          a+=substr(reftrac,'@',substr(ptsyntax,'@',kie[i].abbr))+', ';
          b+='>'+kie[i].abbr+'<:';
          inc(r);
          end;
setlength(a,length(a)-2); a:='    (/ '+a+' /)';  // cut last comma

if (r=0) then out+='! there are no KIEs in selected mechanism, empty arrays are for the compatibility'+_LF;
{$IFDEF USE_PT}
out+='! number of KIE-escorted reactions'+_LF;
out+='  integer, parameter :: NKRREAC = '+inttostr(k)+_LF;
if (r>0) then                                    // no KIE in the system = zero-sized array (for compatibility)
   out+='  integer, parameter :: KRREAC(NKRREAC) = &'+_LF+wrap(a,4,4)
else
    out+='  integer, parameter :: KRREAC(NKRREAC+1) = (/ NKRREAC+1 /)';
out+=_LF;
{$ENDIF}
a:=''; s:=0;
// inspecting r-ns list for tagged species in KIE
for i:=1 to _eqs do
    with eqs[i] do
         if (pos(abbr,b)>0) then   // if in KIE list
            begin
            for j:=1 to 2 do
                if (in_tsl(educ[j])) and (pos('>'+educ[j]+'<',c)=0) then
                   begin
                   a+=substr(tagtrac,'@',educ[j])+', ';
                   c+='>'+educ[j]+'<';
                   inc(s);
                   end;
            for j:=1 to _prod do
                if (in_tsl(prod[j])) and (pos('>'+prod[j]+'<',c)=0) then
                   begin
                   a+=substr(tagtrac,'@',prod[j])+', ';
                   c+='>'+prod[j]+'<';
                   inc(s);
                   end;
            end;

setlength(a,length(a)-2); a:='    (/ '+a+' /)'+_LF;  // cut last comma + in braces

out+='! indices of KIE-related species'+_LF;
out+='  integer, parameter :: NKRSPEC = '+inttostr(s)+_LF;
if (s>0) then                   // no KIE in the system = zero-sized array (for compatibility)
   out+='  integer, parameter :: KRSIND(NKRSPEC) = &'+_LF+wrap(a,4,4)
else
    out+='  integer, parameter :: KRSIND(NKRSPEC+1) = (/ NKRSPEC+1 /)'+_LF;


out+=_LF;
// isotope exchange part

// counting active IEX in the mech
r:=0; s:=0; a:=''; b:=''; c:='';
for i:=1 to _iex do
    if (iex[i].imec) then
       begin
       b+=substr(reftrac,'@',substr(ptsyntax,'@',iex[i].abbr))+', ';
       inc(r);
       end;

for i:=1 to _iesl do
    a+=substr(tagtrac,'@',tsl[iesl[i]].spec)+', ';

setlength(a,length(a)-2); a:='    (/ '+a+' /)';  // cut last comma + in braces
setlength(b,length(b)-2); b:='    (/ '+b+' /)';

if (r=0) then out+='! there is no isotope exchange in selected mechanism, empty arrays are for the compatibility'+_LF;
{$IFDEF USE_PT}
out+='! number of isotope exchange reactions'+_LF;
out+='  integer, parameter :: NIXREAC = '+inttostr(r)+_LF;
if (r>0) then                                    // no KIE in the system = zero-sized array (for compatibility)
   out+='  integer, parameter :: IXREAC(NIXREAC) = &'+_LF+wrap(b,4,4)
else
    out+='  integer, parameter :: IXREAC(NIXREAC+1) = (/ NIXREAC+1 /)';
out+=_LF;
{$ENDIF}
out+='! indices of IEX-related species'+_LF;                                           // '+isoelem+'
out+='  integer, parameter :: NIXSPEC = '+inttostr(_iesl)+_LF;                                      // '+isoelem+'
if (_iesl>0) then                   // no KIE in the system = zero-sized array (for compatibility)
   out+='  integer, parameter :: IXSIND(NIXSPEC) = &'+_LF+wrap(a,4,4)                           // 2x '+isoelem+'
else
    out+='  integer, parameter :: IXSIND(NIXSPEC+1) = (/ NIXSPEC+1 /)'+_LF;                            // 3x '+isoelem+'


imcom_make_kie_iex_relspecs:=out;

end;

// -----------------------------------------------------------------------------
// outputs an array of quantities of isot. element atoms in each species
// prefix: D for doubling, T for tagging
function imcom_make_atom_qties : ansistring;
var i : word;
    a : ansistring;

begin

a:='';
for i:=1 to _utsl do
    a:=a+inttostr(tsl[utsl[i]].qatm)+', ';
setlength(a,length(a)-2); a:='    (/ '+a+' /)';

a:=_LF+
   '! quantity of '+isoelem+' atoms in corresponding tagged species'+_LF+
   '  integer, parameter :: {%NQATOM}({%NSPEC}) = &'+_LF+
   wrap(a,4,4);

imcom_maKe_reps(a);
imcom_make_atom_qties:=a;

end;

// imcom_parse_proc ------------------------------------------------------------
// parses EQN processing code, place to handle KIE
function imcom_parse_proc(fname : string) : ansistring;
var t : text;
    a, u, p, out :  ansistring;
begin

write(_LF+'  imcom_parse_proc(',fname,'): ');

imcom_check_files_exist([fname]);

assign(t,fname);    // processing file
reset(t);

out:='';            // expected output
while not(eof(t)) do
      begin

      readln(t,a);
      imcom_make_reps(a);

      if pos('{>',a)>0 then
         begin
         // get control property and value letter    (e.g. ATOM:C or ISO:2 or SPEC:Ch3Br etc.)
         p:=upcase(trim(imcom_ext4marks(a,'{>',':')));
         u:=trim(imcom_ext4marks(a,':','}'));
//       out+='! ----- imcom ----- condition in effect ('+p+':'+u+'): start -----'+_LF;
         if not(imcom_condition(p,u)) then
            begin
            repeat
              readln(t,a);
//            out+='! '+a;
            until ( ((pos(('{<'),a)>0) and (pos(uppercase(p),uppercase(a))>0)) or (eof(t)) );
//          out+='! ---- imcom ----- condition in effect ('+p+':'+u+'): stop -----'+_LF;
            end;
         end
      else
          if pos('{<',a)>0 then
             begin
             {do nothing}
//           out+='! ---- imcom ----- condition in effect ('+p+':'+u+'): stop -----'+_LF;
             end
          else
              out:=out+a+_LF;
      end;

close(t);

writeln('done');

imcom_parse_proc:=out;

end;


// imcom_parse_eq_str_arrlen ---------------------------------------------------
// parses given string and expands with spaces (or shrinks) any string
// expression occurence (i.e. 'I12DOsa') to the length of the type nstr
// (made foe fortran compilers compatibility
function imcom_parse_eq_strarr_len(s : ansistring) : ansistring;
var o : ansistring;
    i, l, z : longint;
    m : boolean;

begin

z:=sizeof(nstr)-1;
o:='';
m:=false; l:=0;
for i:=1 to length(s) do
    begin

    if s[i]='''' then
       begin
       m:=not(m);       // m is a flag "inside of the string"
       if m then l:=0   // l is a length of the element
            else
                while (l<=z) do
                      begin o+=' '; inc(l); end;

       end;

    if (l<=z) or (s[i]='''') or not(m) then
       o+=s[i];
    if (m) then inc(l);

    end;

imcom_parse_eq_strarr_len:=o;

end;


// --- update_ptracs -----------------------------------------------------------
// adds tagged reactions abbrs. to the inter-conf. tagged reacs. list
procedure imcom_update_ptracs;

{$IFDEF USE_PT}
procedure add_ptrac(abbr : string);
begin
if (pos(':'+abbr+'>',ptracs_intr)=0) then
   begin
   ptracs_intr+=':'+abbr+'>';           // adding to the PTs list
   inc(_ptracs_intr);
   end;
if (pos(':'+abbr+'>',ptracs_conf)=0) then
   begin
   ptracs_conf+=':'+abbr+'>';           // adding to the PTs list
   inc(_ptracs_conf);
   end;
end;

var i, j, k : integer;
    has_pt : boolean;                   // flag showing if the eqn has already a PT
begin

for i:=1 to _eqs do
    if (eqs[i].itag) then               // if tagged, should have a passive tracer added
       begin
       has_pt:=false;
       for j:=1 to eqs[i]._prod do
           if (eqs[i].prod[j]=substr('@',eqs[i].abbr,ptsyntax)) then
              begin has_pt:=true; break; end;
       if not(has_pt) then
          add_ptrac(eqs[i].abbr);
       end;

// additional tracers for KIE monitoring
{$IFDEF USE_PT_KIE}
for k:=1 to _kie do
    if (kie[k].imec) then       // if KIE exist in this configuration
       with kie[k] do
            for j:=1 to _isos do
                if not( in_tsl(eqs[kie[k].eqno].educ[1]) and in_tsl(eqs[kie[k].eqno].educ[2]) ) then    // not a quadrupl. reaction
                   add_ptrac(abbr+cfgname+clsname[j])
                else
 {$IFNDEF TAG_EXPL}
                    begin            // in case both educts are tagged (quadrupled equation = quad. kie PTs)
                    add_ptrac(kie[k].abbr+cfgname+clsname[j]+'i1');
                    add_ptrac(kie[k].abbr+cfgname+clsname[j]+'i2');
                    end;
 {$ELSE}
                    for i:=1 to _isos do  // in case of explicit doubling (_iso-replicated equations with KIE)
                        add_ptrac(kie[k].abbr+cfgname+clsname[j]+'e'+inttostr(i));
 {$ENDIF}
{$ENDIF}

writeln('imcom_update_ptracs: done [',_ptracs_conf,' configuration, ',_ptracs_intr,' total PTs]');
{$ELSE}  // USE_PT
begin
{$ENDIF} // USE_PT
end;


// imcom_ptracs_list (returns string with the passive tracers list -------------
// outputs the passive tracers list according to specimen with @, $
function imcom_make_ptracs_list(ptstr, sample : string) : ansistring;
var tmp, smp, out : ansistring;
    no : integer;
begin
out:='';
tmp:=ptstr;
no:=1;
while (length(tmp)>0) do
      begin
      // one ptrac entry is like :G4410>
      smp:=copy2symbdel(tmp,':');   // rem ':'
      smp:=copy2symbdel(tmp,'>');   // copy to '>'

      smp:=substr(sample,'@',smp);
      smp:=substr(smp,'$',inttostr(no));

      out+=smp;
      inc(no);
      end;
out:=trim(out);                // cut leading spaces
imcom_make_ptracs_list:=out;
end;

// imcom_make_ptracs (returns string with the passive tracers indices defs) ----
// in: boolean indicating if the list is created for the inter. former
function imcom_make_ptracs(is_intr : boolean) : ansistring;
var a, out, pts : ansistring;
    m : nstr;
{$IFDEF USE_PT}
    i, j, _pts : integer;
const l : char = 'T';
{$ENDIF}
begin

out:='';

{$IFDEF USE_PT}
if (is_intr) then
   begin _pts:=_ptracs_intr; pts:=ptracs_intr; m:=''; end
else
    begin _pts:=_ptracs_conf; pts:=ptracs_conf; m:=isoelem; end;

out+='! number of tagged reactions'+_LF;
out+='  integer, parameter :: N'+l+m+'REAC = '+inttostr(_pts)+_LF;
out+='! added passive tracers indices'+_LF;
out+='  integer, parameter :: '+l+'RPT'+m+'IND(N'+l+m+'REAC) = &'+_LF;
out+='    (/ ';
a:=imcom_make_ptracs_list(pts,' '+substr(sisyntax,'@',ptsyntax)+',');
setlength(a,length(a)-1); // cut last comma
out+=wrap(a,7,7);
out+=' /)'+_LF;

// for the interconfiguration former the internal reaction numbers can be added
if (is_intr) then
   begin
   j:=0; a:='';
   for i:=1 to _eqs do
       if (eqs[i].etag>0) then
          a+=inttostr(i)+', ';

   setlength(a,length(a)-2); // cut last comma

   out+=_LF+'! tagged reactions internal nos.'+_LF;
// hopefully, N#RININD equals to N#REAC, otherwise diagnose internal problem!
// out+='  integer, parameter :: N'+l+'RININD = '+inttostr(j)+_LF;
// out+='  integer, parameter :: '+l+'RININD(N'+l+'RININD) = &'+_LF;
   out+='  integer, parameter :: '+l+'RININD(N'+l+m+'REAC) = &'+_LF;
   out+='    (/ '+wrap(a,7,7)+' /)'+_LF;
   end;

{$ENDIF}

imcom_make_ptracs:=out;
end;

// *** OUTDATED, TO BE REMOVED ***
// imcom_resetPTs --------------------------------------------------------------
// makes PTs (production tracers) initialization (reset to zero) code
// using specimens for index-name and value-expression (like ind_@ and C(ind_@)
function imcom_make_resetPTs(valpref, indpref : string) : ansistring;
var i : word;
    a : ansistring;
{$IFDEF USE_PT_KIE}
    b : ansistring;
    j, e : word;
 {$IFDEF TAG_EXPL}
    x : word;
 {$ENDIF}
{$ENDIF}

begin
{$IFDEF USE_PT}

if indpref='' then indpref:='ind_'+ptsyntax;
if valpref='' then valpref:='C('+indpref+')';

a:='';

// unset USE_IF_IN_RESETPT to optimize (attention!)

a:=a+'! - [imtag/double] generated ------------------------------------------------'+_LF+_LF;
a:=a+'! production tracers initialization ('+indpref+', '+valpref+')'+_LF+_LF;

{$IFNDEF TAG_EXPL}
for i:=1 to _eqs do
    if (eqs[i].itag) then
       begin
       a:=a+'    ';
{$IFDEF USE_IF_IN_RESETPT}
       a:=a+'IF ('+substr(indpref,'@',eqs[i].abbr)+' /= 0) ';
{$ENDIF}
       a:=a+substr(valpref,'@',eqs[i].abbr)+' = 0.0_dp'+_LF;
       end;
{$ENDIF}

{$IFDEF USE_PT_KIE}        // if monitoring KIE, all isotopologues (only for doubling)
a:=a+_LF+'  ! reactions with KIE'+_LF;
b:='';                         // added reactions list
for i:=1 to _kie do
    for e:=1 to _eqs do
        if (eqs[e].itag) then
           if (eqs[e].abbr=kie[i].abbr) and (pos('>'+kie[i].abbr+'<',b)=0) then
              begin
              for j:=1 to _isos do
                  if not(in_tsl(eqs[e].educ[1]) and in_tsl(eqs[e].educ[2])) then
                     begin
                     a:=a+'    ';
 {$IFDEF USE_IF_IN_RESETPT}
                     a:=a+'IF ('+substr(indpref,'@',kie[i].abbr+clsname[j]{+kie[i].isot}+' /= 0) ';
 {$ENDIF}
                     a:=a+substr(valpref,'@',kie[i].abbr+clsname[j]{+kie[i].isot})+' = 0.0_dp'+_LF;
                     end
                  else
 {$IFNDEF TAG_EXPL}
                      begin            // in case both educts are tagged (quadrupled equation = quad. kie PTs)
                      a:=a+'    ';
  {$IFDEF USE_IF_IN_RESETPT}
                      a:=a+'IF ('+substr(indpref,'@',kie[i].abbr+clsname[j]+'i1')+' /= 0) ';
  {$ENDIF}
                      a:=a+substr(valpref,'@',kie[i].abbr+clsname[j]+'i1')+' = 0.0_dp'+_LF;

                      a:=a+'    ';
  {$IFDEF USE_IF_IN_RESETPT}
                      a:=a+'IF ('+substr(indpref,'@',kie[i].abbr+clsname[j]+'i2')+' /= 0) ';
  {$ENDIF}
                      a:=a+substr(valpref,'@',kie[i].abbr+clsname[j]+'i2')+' = 0.0_dp'+_LF;
                      end;
 {$ELSE}
                      for x:=1 to _isos do  // in case of explicit doubling (_iso-replicated equations with KIE)
                          begin
                          a:=a+'    ';
  {$IFDEF USE_IF_IN_RESETPT}
                          a:=a+'IF ('substr(indpref,'@',kie[i].abbr+clsname[j]+'e'+inttostr(x))+' /= 0) ';
  {$ENDIF}
                          a:=a+substr(valpref,'@',kie[i].abbr+clsname[j]+'e'+inttostr(x))+' = 0.0_dp'+_LF;
                          end;
 {$ENDIF}
              b:=b+'>'+kie[i].abbr+'<';
              end;
{$ENDIF}

imcom_make_resetPTs:=a;
{$ENDIF}
end;


// imcom_make_tracdecl ---------------------------------------------------------
// outputs tracers declaration code
function imcom_make_tracdecl(tagtrac : string) : ansistring;

var i, j, s : word;
    a, b, out : ansistring;

begin

// default "tagged" tracers naming
if tagtrac='' then tagtrac:='ind_t@';

out:='! data declaration'+_LF;
out+=_LF;
out+='! number of tagged species'+_LF+
     '  integer, parameter :: {%NSPEC} = '+inttostr(_utsl)+_LF+
     '! number of classes (isotopologues)'+_LF+
     '  integer, parameter :: {%NCLASS} = '+inttostr(_isos)+_LF+
     '  integer, parameter :: {%NISO} = '+inttostr(_isos)+_LF+
     '! number of "fixed" specs'+_LF+
     '  integer, parameter :: {%NFIX} = '+inttostr(_fixs)+_LF;

out+=_LF;
out+='! isotopologues tracers indices'+_LF;
for s:=1 to _utsl do
    with tsl[utsl[s]] do
         begin
         out+='  integer, parameter :: '+substr(tagtrac,'@',spec)+' = '+inttostr(s);
         a:=''; for j:=1 to _isos do 
           if (ifix[j]) then a+=' '+clsname[j];
         if (a<>'') then out+='          ! classes'+a+'indicated as "fixed"';
         out+=_LF;
         end;
out+=_LF;

{$IFDEF ADD_DUMIND}
out+='! dummy indices of species present TSL but not in the mechanism (resolves dependency problems)'+_LF;
for s:=1 to _tsl do
    if (not(in_utsl(tsl[s].spec))) then
       out+='  integer, parameter :: '+substr(tagtrac,'@',tsl[s].spec)+' = 0'+_LF;
out+=_LF;
{$ENDIF}


// --- species indices ---
out+='! corresponding to regular(0) <-> tagged (1, 2, ...) species indices '+_LF;
out+='  integer, parameter :: {%RSIND}({%NSPEC},0:{%NCLASS}) = reshape( &'+_LF;

a:='    (/ ';

// regulair
for i:=1 to _utsl do
    a+=substr(sisyntax,'@',tsl[utsl[i]].spec)+', ';
// en toute, esthetique
a+='&';

// isotopologues
for j:=1 to _isos do
    begin
    for i:=1 to _utsl do
        a+=substr(sisyntax,'@',tsl[utsl[i]].isos[j])+', ';

    // en toute, esthetique
    a+='&';
    end;
setlength(a,length(a)-3);  // cut last comma and ampersand
a+=' /), (/ {%NSPEC}, 1+{%NCLASS} /) )';

out+=wrap(a,7,7)+_LF;
out+=_LF;


{$IFDEF USE_PT_UPL}
// --- unacccounted production/loss ---
out+='! total unaccounted production/loss PTs indices'+_LF;
out+='  integer, parameter :: {%UPLIND}(2,0:{%NCLASS}) = reshape( &'+_LF;
// regulair
a:='    (/ '+substr(sisyntax,'@',substr(ptpsyntax,'@','U'+cfgname+isoelem))+', '+substr(sisyntax,'@',substr(ptlsyntax,'@','U'+cfgname+isoelem));
// isotopiques
for j:=1 to _isos do
    a+=', '+substr(sisyntax,'@',substr(ptpsyntax,'@','U'+cfgname+clsname[j]))+', '+substr(sisyntax,'@',substr(ptlsyntax,'@','U'+cfgname+clsname[j]));
// en toute, esthetique
a+=' /), (/ 2, 1+{%NCLASS} /) )';
out+=wrap(a,7,7);

out+=_LF;
{$ENDIF}

// --- "fixed" species list ---
a:='';
out+=_LF;
out+='! list of "fixed" species indices'+_LF;
if (_fixs=0) then
   begin
   out+='! there is no "fixed" species in selected mechanism, zero array is for the code compatibility'+_LF;
// out+='  integer, parameter :: {%FSIND}({%NFIX}+1,0:{%NCLASS}) = 0'+_LF; // as a parameter yields boundary check error with GCC :(
   out+='  integer            :: {%FSIND}({%NFIX}+1,0:{%NCLASS}) = 0'+_LF;
   end
else
    begin
    out+='  integer, parameter :: {%FSIND}({%NFIX},0:{%NCLASS}) = reshape( &'+_LF;
      a:='    (/ ';
    for s:=1 to _utsl do
        begin
        i:=0; b:='';
        for j:=1 to _isos do
            if (tsl[utsl[s]].ifix[j]) then 
               begin
               i:=i+1;
               b+=substr(sisyntax,'@',tsl[utsl[s]].isos[j])+', ';
               end
            else
                b+='0, ';
        if (i=0) then continue; // skip species without fixed classes
        a+=substr(tagtrac,'@',tsl[utsl[s]].spec)+', '+b;
        end;
    if (_fixs>0) then setlength(a,length(a)-2);  // cut last comma, if there's any
    a+=' /), (/ {%NFIX}, 1+{%NCLASS} /) )';
    end;

out+=wrap(a,7,7);

out+=_LF+imcom_make_atom_qties+_LF;

{$IFDEF USE_PT}
out+=_LF+imcom_make_ptracs(false)+_LF;
{$ENDIF}

{$IFDEF USE_KRSIND}
out+=_LF+imcom_make_kie_iex_relspecs(sisyntax,tagtrac)+_LF;
{$ENDIF}

imcom_make_reps(out);
imcom_make_tracdecl:=out;

end;


// imcom_make_specslist ---------------------------------------------------------
// outputs species list from according to a sample
function imcom_make_specslist(sample : string) : ansistring;

var i : integer;
    a, b, out : ansistring;

begin
out:='';

b:=sample;
imcom_make_reps(b);

for i:=1 to _utsl do
    begin
    a:=substr(b,'$',inttostr(i));
    a:=substr(a,'@',tsl[utsl[i]].spec);
    out+=a;
    end;

imcom_make_specslist:=out;
end;



// imcom_x0 ----------------------------------------------------
// code for x0 initialization procedure
function imcom_make_x0(classno, expstr : string) : ansistring;
var i, c, cmin, cmax : integer;
    out, a : ansistring;
begin

// getting the class # (Nth - for Nth minor isotopologue or class)
if (pos('#',classno)>0) then
   begin cmin:=1; cmax:=_isos; end
else
    begin
    cmin:=strtointdef(trim(classno),1);         // default = 1
    cmax:=cmin;
    end;

out:='  ! TODO: specify the initial values here:'+_LF;
out+=_LF;

for c:=cmin to cmax do
    begin

    if (c>_isos) then
       begin
       write('<!> ');
       imcom_make_x0:='  ! imcom_make_x0('+trim(classno)+','+expstr+'): minor class # exceeds number of given in this configuration <may be normal>';
       exit;
       end;

    if (c=1) and (_isos>1) and (itransm=1) then
       begin
       write('<!> ');
       out+='  ! imcom_make_x0('+trim(classno)+','+expstr+'): requested class #1 for configuration with '+inttostr(_isos)+' classes <may be normal>'+_LF+_LF;
       end;

    // by default: MECCA indices
    if (expstr='') then expstr:='INIT_'+classno+'(ind_@) = $'; // something to have by default

    out+='  ! class: '+clsname[c]+_LF;
    for i:=1 to _utsl do
        if (tsl[utsl[i]].ival[c]<>'') then
           out+=substr( substr( substr(expstr,'@',tsl[utsl[i]].spec),
                                              '$',tsl[utsl[i]].ival[c]),
                                              '#',inttostr(c))+_LF;

    out+=_LF;
    end;

a:='';
for i:=1 to _tsl do          // assuming #1 - 12C, #2 - 13C, tagging 12C/13C
    if not(is_usedspec(tsl[i].spec)) then
       a+='! '+tsl[i].spec+', ';
setlength(a,min(length(a)-2,0));
if (length(a)>0) then
  out+=_LF+'! following species were not found in current selected mechanism:'+_LF+wrap(a,0,2);

imcom_make_x0:=out;

end;



// imcom_f0 ----------------------------------------------------
// code for f0 adjustment procedure - for fixed specs only
function imcom_make_f0(classno, expstr : string) : ansistring;
var i, c, cmin, cmax : integer;
    out, a : ansistring;
begin

// getting the class # (Nth - for Nth minor isotopologue or class)
if (pos('#',classno)>0) then
   begin cmin:=1; cmax:=_isos; end
else
    begin
    cmin:=strtointdef(trim(classno),1);         // default = 1
    cmax:=cmin;
    end;

out:='  ! TODO: specify the initial values here:'+_LF;
out+=_LF;

for c:=cmin to cmax do
    begin

    if (c>_isos) then
       begin
       write('<!> ');
       imcom_make_f0:='  ! imcom_make_f0('+trim(classno)+','+expstr+'): minor class # exceeds number of given in this configuration <may be normal>';
       exit;
       end;

    if (c=1) and (_isos>1) then
       begin
       write('<!> ');
       out+='  ! imcom_make_f0('+trim(classno)+','+expstr+'): requested class #1 for configuration with '+inttostr(_isos)+' classes <may be normal>'+_LF+_LF;
       end;

    // by default: MECCA indices
    if (expstr='') then expstr:='INIT_'+classno+'(ind_@) = $'; // something to have by default

    out+='  ! class: '+clsname[c]+_LF;
    for i:=1 to _utsl do
        if (tsl[utsl[i]].ival[c]<>'') and (tsl[utsl[i]].ifix[c]) then
           out+=substr( substr( substr(expstr,'@',tsl[utsl[i]].spec),
                                              '$',tsl[utsl[i]].ival[c]),
                                              '#',inttostr(c))+_LF;

    out+=_LF;
    end;

a:='';
for i:=1 to _tsl do          // assuming #1 - 12C, #2 - 13C, tagging 12C/13C
    if not(is_usedspec(tsl[i].spec)) then
       a+='! '+tsl[i].spec+', ';
setlength(a,min(length(a)-2,0));
if (length(a)>0) then
  out+=_LF+'! following species were not found in current selected mechanism:'+_LF+wrap(a,0,2);

imcom_make_f0:=out;

end;

// - util ----------------------------------------------------------------------

// -- checks and corrects the duplicate equations
procedure imcom_check_eqn_dupes(fname, dumspec : string);

function is_dupe_eq(eq1, eq2 : ansistring) : boolean;
var sl1, sl2 : tstringlist;
    sample : ansistring;
    i : integer;
const _esd_ = ['+','-'];        // delimiter of the members of the side of eq

begin

sl1:=tstringlist.create;
sl2:=tstringlist.create;

// processing left side
sl1.clear; sl2.clear;

sample:=copy(eq1,1,pos('=',eq1)-1);
for i:=1 to wordcount(sample,_esd_) do
    sl1.add( extractword(i,sample,_esd_) );

sample:=copy(eq2,1,pos('=',eq2)-1);
for i:=1 to wordcount(sample,_esd_) do
    sl2.add( extractword(i,sample,_esd_) );

sl1.sort; sl2.sort;

result:=sl2.equals(sl1);

if (result) then
   begin

   // processing right side
   sl1.clear; sl2.clear;
   sample:=copy(eq1,pos('=',eq1)+1,length(eq1)-pos('=',eq1));
   for i:=1 to wordcount(sample,_esd_) do
       sl1.add( extractword(i,sample,_esd_) );

   sample:=copy(eq2,pos('=',eq2)+1,length(eq2)-pos('=',eq2));
   for i:=1 to wordcount(sample,_esd_) do
       sl2.add( extractword(i,sample,_esd_) );

   sl1.sort; sl2.sort;

   result:=result and sl1.equals(sl2);
   end;

//if (result) then
//   begin
//   writeln('-dupe----------');
//   writeln(eq1);
//   writeln(eq2);
//   end;

sl1.destroy;
sl2.destroy;

end;

var f, o : text;
    e : array[1..max_eqs*max_isos] of record
        eq : ansistring;
        line, corr : integer;
        end;
    a : ansistring;
    l, p, i, j : integer;

begin

write('check_eqn_dupes(',fname,'): ');

assign(f,fname);
reset(f);

p:=0; l:=0;
while not(eof(f)) do
      begin
      inc(l);
      readln(f,a);
      if imcom_is_eqn(a) then
         begin
         inc(p);
         e[p].eq:=imcom_ext4marks(a,'>',':');          // copier l'equation
         e[p].eq:=uppercase(delspace(trim(e[p].eq)));  // >JE+SUIS+UNE=EQUATION+0.5MAINTENANT:
         e[p].line:=l;
         e[p].corr:=0;
         end;
      end;
close(f);

// bubblesort comparison & correction
l:=1;
for i:=1 to p-1 do
    begin

    // most time-taking part, showing runner
    show_runner(p div 32,i+j);

    for j:=i+1 to p do
        if is_dupe_eq(e[j].eq,e[i].eq) then
           begin
           inc(e[j].corr);
           e[j].eq+=uppercase('+'+dumspec);
           end;
    end;

// apporter les modifications
assign(o,fname+'.$$$'); rewrite(o);
assign(f,fname); reset(f);

j:=1; // indice a e[]
l:=0; // indice a f
p:=0;
while not(eof(f)) do
      begin
      inc(l);
      readln(f,a);
      if (e[j].line=l) then
         begin

         if (e[j].corr>0) then
            begin
            insert('+ '+inttostr(e[j].corr)+' '+dumspec+' ',a,pos(':',a));
            inc(p);
            end;

         inc(j);
         end;
      writeln(o,a);
      end;
close(o); close(f);
erase(f); rename(o,fname);

if (p>0) then
   write('corrected ',p,' identical equations')
else
    write('no identical equations found');

writeln;

end;

// -- returns the parameters-defines block
function imcom_parameters : ansistring;
var s : ansistring;
    i : integer;
const ds : ansistring = '#define ';
begin

s:='! configuration parameters'+_LF+_LF;

// tagging name
s+=ds+'TAG'+_LF;
s+=ds+'tag_'+cfgname+_LF;
 // explicit/implicit doubling
 {$IFDEF TAG_EXPL}
 s+=ds+'TAG_EXPLICIT';
 {$ELSE}
 s+=ds+'TAG_IMPLICIT';
 {$ENDIF}
s+=_LF;

s+=_LF;

// abbr
s+=ds+'CONF_'+cfgname+_LF;

// atom
s+=ds;
if (isoelem<>'') then s+='ATOM_'+isoelem else s+='NO_ATOM';
s+=_LF;

// classes no
s+=ds+'CLASSES_'+inttostr(_isos)+_LF;

// classes names
for i:=1 to _isos do
    s+=ds+'CLASS_'+clsname[i]+_LF;

s+=_LF;

// kie presence
s+=ds;
if (_kie>0) then s+='KIE' else s+='NO_KIE';
s+=_LF;

// source specification presence
s+=ds;
if (_src>0) then
   s+='SRC' else s+='NO_SRC';
s+=_LF;

// isotope exchange reactions presence
s+=ds;
if (_iex>0) then
   s+='IEX' else s+='NO_IEX';
s+=_LF;

s+=_LF;

// usage of passive tracers flag
{$IFDEF USE_PT}
s+=ds+'USE_PT'+_LF;
 {$IFDEF USE_PT_KIE}
 s+=ds+'USE_PT_KIE'+_LF;
 {$ENDIF}
{$ENDIF}

// usage of reaction numbers flag
{$IFDEF USE_RN}
s+=ds+'USE_RN'+_LF;
{$ENDIF}

if (cparams<>'') then
   begin
   i:=1;
   s+=_LF+'! user configuration parameters'+_LF;
   while (extractword(i,cparams,_delims)<>'') do
         begin
         s+=_LF+ds+extractword(i,cparams,_delims);
         inc(i);
         end;
   s+=_LF;
   end;

// replacements made
s+=_LF+'! replacements:'+_LF+'!';
for i:=1 to _imcom_reps do
    s+=_LF+'! '+imcom_reps[i,1]+' is '+imcom_reps[i,2];

imcom_parameters:=s;

end;

// - inter-configuration part --------------------------------------------------

// maximum number of different configurations allowed
const max_conf = 5;

var conf : array[1..max_conf] of string;  // configuration-ID
   _conf : integer;                       // no. of configurations

// imcom_make_configslist ------------------------------------------------------
// make a list of samples with configurations names (for inter-conf. modules)
function imcom_make_configslist(sample : ansistring) : ansistring;
var i : integer;
    out : ansistring;
begin
out:='';
for i:=1 to _conf do
    out+=substr(sample,'@',conf[i])+_LF;
setlength(out,length(out)-length(_LF));
imcom_make_configslist:=out;
end;


// -----------------------------------------------------------------------------
// imcom_init: initialization of some inter-configuration parameters (to safe 0s)
procedure imcom_init;
begin
_conf:=0;
_form_intr:=0;
{$IFDEF USE_PT}
   _ptracs_intr:=0; ptracs_intr:='';
{$ENDIF}
// default time output format
defaultformatsettings.shortdateformat:='DD-MMM-YYYY';
defaultformatsettings.longtimeformat:='hh:nn:ss';
// runner counter
runner_cnt:=0;
runner_chr:=1;
end;



// EOF
