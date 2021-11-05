
// = imtag.pas =======================================================
// mz_sg_20200715 (gromov@tonnerre.mpic.de)
//
// (isotopic) mechanism tagging tool, version 2.6
//
// im-tag main source module
// requirements: free pascal compiler (fpc) ver. 3 or higher
//
// [S.Gromov, MPIC, 2007-2020]
// ===================================================================

program imtag;

// - conditional defines -------------------------------------------

{$DEFINE xTAG_EXPL}         // explicit tagging - original reactions are replaced
                            // can be defined during compilation with -dTAG_EXPL switch

{$DEFINE xUSE_PT}           // using passive tracers, may be OFF
{$DEFINE xUSE_PT_KIE}       // defines whether to add PTs to all KIE-reactions (for monitoring)
// {$DEFINE  USE_PT_UPL}    // defines whether to add PTs for unaccounted elemental/molecular comp. production/loss  <-- set via Makefile
{$DEFINE  USE_KRSIND}       // defines whether to output the list of KIE-related specs indices (for correction)

{$DEFINE xPRECISE_FRACISO_IEX} // defines whether to treat IEX in case of "fractional" isotopologue
                               // (i.e. 1 class) very carefully (i.e. weighting of reaction rate to regular-rare)

{$DEFINE  USE_DKRATE}       // defines whether to put reaction rates in the eqn header as a
                            // separate variables and use them in multiple tagged reactions

{$DEFINE  IEX_ONWARD_RATES} // defines whether isotope exchange rates are given as onward for the
                            // first reactant, or referred to the regular reaction (affects the rate coefficient)

{$DEFINE xIEX_REGREF}       // defines whether to enable creation of regular reactions as a reference
                            // for the isotope exchange reactions

// -----------------------------------------------------------------

{$I imcom.inc}              //  linking common utils/types include file for im-tag

// -----------------------------------------------------------------

{$IFDEF USE_DKRATE}
const prev_ptracs_intr : ansistring = ''; // previous value of ptracs_intr
      dkratesyntax = 'tag_k@';            // naming of vars for duped r.rate coeffs
{$ENDIF}

const dsrcsubsyntax = 'tag_src_f@';       // naming of vars for source substit. fracs

const stoiformat = '0.#######';
      dummyspc = 'Dummy';

// main tagging routine
procedure tagit(fname : string);

function giveiso(spec : string; ni : integer; abbr : nstr = '') : string;
var ns, s, h : integer;
begin

h:=0;
ns:=no_tsl(spec);      // get a number of a species in tsl
if (ns=0) or (ni<1) or (ni>_isos) then
   giveiso:=spec
else
    if (abbr<>'') then // if abbreviation is given, we search/apply shifts
       for s:=1 to _shf do
           if (pattmatch(abbr,shf[s].eqn,true) and pattmatch(spec,shf[s].spc,true)) then
              begin
              case shf[s].mode of
                   shf_next_acc  : ni:=min(ni+1,_isos);
                   shf_prev_acc  : ni:=max(ni-1,1);
                   shf_next_loss : inc(ni);
                   shf_prev_loss : dec(ni);
              end; {case}
              if (ni<1) or (ni>_isos) then  // in case of loss
                 giveiso:=dummyspc
              else
                  giveiso:=tsl[ns].isos[ni];
              //writeln(' <*> caught ',spec,' in r-n ',abbr,' mode: ',shf[s].mode);
              exit; //inc(h);
              end;

if (h=0) then
   giveiso:=tsl[ns].isos[ni];

if (h>1) then
   writeln(' <!> warning: giveiso: multiple (',h,') shift conditions met for ',spec,
           ' in r-n ',abbr,' (you should check tagging cfg, [SHF] section');

end;

var f : textfile;     // output file

// creation of inline-part for updating of a
// rate constants coefficients for a specified sources
procedure addinline4subs;

var i, j, k : word;
    a : ansistring;
    srcspecs : array[1..max_tsl] of nstr;
   _srcspecs : integer;
    gotcha : boolean;
begin

// making the list of all "source" species used in source specification
_srcspecs:=0;
fillchar(srcspecs,sizeof(srcspecs),0);
for i:=1 to _eqs do
    if ((eqs[i].itag) and (eqs[i].isrc)) then
    for j:=1 to src[eqs[i].nsrc]._trans do
        begin
        gotcha:=false;              // checking if a spec in the list already
        for k:=1 to _srcspecs do
            if (src[eqs[i].nsrc].trans[j].src=srcspecs[k]) then gotcha:=true;
        if (not(gotcha) and is_usedspec(src[eqs[i].nsrc].trans[j].src)) and
           ((src[eqs[i].nsrc].trans[j].src<>eqs[i].educ[1]) and   // if the source is on the lef side of eq, no fraction calc is needed
            (src[eqs[i].nsrc].trans[j].src<>eqs[i].educ[2])) then
           begin
           inc(_srcspecs);
           srcspecs[_srcspecs]:=src[eqs[i].nsrc].trans[j].src;
           end;
        end;

{$IFDEF PRECISE_FRACISO_IEX}
// source species for IEX for "fractional" isotopologues
if not(_isos>1) then
   for i:=1 to _eqs do
       if ((eqs[i].itag) and (eqs[i].iiex)) then
         with iex[eqs[i].niex] do
          for j:=1 to 2 do
              begin
              gotcha:=false;              // checking if a spec in the list already
              for k:=1 to _srcspecs do
                  if (tsl[exspec[j]].spec=srcspecs[k]) then gotcha:=true;
              if (not(gotcha) and is_usedspec(tsl[exspec[j]].spec)) then
                 begin
                 inc(_srcspecs);
                 srcspecs[_srcspecs]:=tsl[exspec[j]].spec;
                 end;
              end;
{$ENDIF}

// in case we've nothing to add
if (_srcspecs=0) then exit;

// GLOBAL part
writeln(f,'#INLINE F90_GLOBAL');
writeln(f,'! ---------------- [',tagname,'] - inline code ----------------');
writeln(f,'! reaction rates modifiers according to the specified source isotopic composition');

writeln(f,'  real(dp) :: '+tagname+'_src_temp         ! variable used for fractions calculation');

for i:=1 to _srcspecs do
    begin
//       writeln(f,'! '+src[i].spec+' isotopologue fractions');
    a:='  real(dp) :: ';
    for j:=1 to _isos do
        a:=a+substr(dsrcsubsyntax,'@',clsname[j]+srcspecs[i])+', ';
    setlength(a,length(a)-2);
    writeln(f,wrap(a,2,4));
    end;
writeln(f,'#ENDINLINE');

// RCONST part
writeln(f,'#INLINE F90_RCONST');
writeln(f,'! ---------------- [',tagname,'] - inline code ----------------');
writeln(f,'! calculation of the isotopologue fractions');
//writeln(f,'!');

for i:=1 to _srcspecs do
    begin
    // total composition
    a:='(';
    if (_isos>1) then
       for j:=1 to _isos do
           a+='C('+substr(sisyntax,'@',clsname[j]+srcspecs[i])+')+'            // substr(sisyntax,'@',clsname[j]+srcspecs[i]) should be replaced with isos[]
    else
        a+='C(ind_'+srcspecs[i]+')+';
    setlength(a,length(a)-1);
    a+=')';
    if not(_isos>1) then
       a+='       ! <!> warn: 1 class used, thus weighting to regular ';
    writeln(f,'  '+tagname+'_src_temp = ',a);
    writeln(f,'  if ('+tagname+'_src_temp .gt. 0.0_dp) then');
    for j:=1 to _isos do
        writeln(f,'    '+substr(dsrcsubsyntax,'@',clsname[j]+srcspecs[i])+' = ',
                    'C(ind_'+clsname[j]+srcspecs[i]+') / '+tagname+'_src_temp');
    writeln(f,'  else');
    for j:=1 to _isos do
        writeln(f,'    '+substr(dsrcsubsyntax,'@',clsname[j]+srcspecs[i])+' = 0.0_dp');
    writeln(f,'  endif');
//  writeln(f,'!');
    end;
writeln(f,'#ENDINLINE');

end;

// creation of inline-part for updating of a
// rate constants coefficients for a tagged equations
procedure addinline4tagrates;
var i : word;
    a : ansistring;
    any2rate : boolean;

begin
{$IFDEF USE_DKRATE}
// detecting if there are reactions not yet processed
any2rate:=false;
for i:=1 to _eqs do
    if (eqs[i].itag) and (eqs[i].etag<2) then
       any2rate:=true;
if not(any2rate) then exit;

// GLOBAL part
writeln(f,'#INLINE F90_GLOBAL');
writeln(f,'! ---------------- [',tagname,'] - inline code ----------------');
writeln(f,'!');
writeln(f,'! reaction rates variables (to be used in tagged reactions)');
a:='  real(dp) :: ';
for i:=1 to _eqs do
    if (eqs[i].itag) and (eqs[i].etag<2) then
//       if (pos(':'+eqs[i].abbr+'>',prev_ptracs_intr)=0) then    // var. might have been added in the previous conf.
          a+=substr(dkratesyntax,'@',eqs[i].abbr)+', ';
setlength(a,length(a)-2);
writeln(f,wrap(a,2,4));
writeln(f,'#ENDINLINE');

// RCONST part
writeln(f,'#INLINE F90_RCONST');
writeln(f,'! ---------------- [',tagname,'] - inline code ----------------');
writeln(f,'! reaction rates variables calculation (to be calculated once per rates evaluation)');
writeln(f,'!');
for i:=1 to _eqs do
    if (eqs[i].itag) and (eqs[i].etag<2) then
//     if (pos(':'+eqs[i].abbr+'>',prev_ptracs_intr)=0) then    // var. might have been added in the previous conf.
          begin
//        writeln(f,'  '+substr(dkratesyntax,'@',eqs[i].abbr)+' = '+trim(imcom_ext4marks(eqs[i].phys,'}',';')+';'));
          // cutting Domenico- and Rolf- style comments, i.e. "{10^ (+-3)}" or "// something" in the r-n rate :-]
            writeln(f,'  '+substr(dkratesyntax,'@',eqs[i].abbr)+' = '+trim(substr(imcom_rem4marks(eqs[i].phys,'{','}'),'//','! ')));
//crashes/  writeln(f,' ! ',trim(delspace1(imcom_ext4marks(eqnfile[eqs[i].eqnf].line,'>',':')))); // r-n ref. helper
          end;
writeln(f,'#ENDINLINE');
{$ENDIF}

{$IFDEF USE_PT}
prev_ptracs_intr:=ptracs_intr;
{$ENDIF}
end;


// function for formatting products
function format_prod(prod : nstr; stoi : real) : string;
var o : string;
begin
if (abs(stoi)<eps_zero) then
   format_prod:=''  // products with null stoi are not output
else
    begin
    o:='';
    if (abs(stoi-1)>eps_zero) then
       o+=formatfloat(stoiformat,stoi)+' '; // stoichiometric coeff. (if not unity)
    o+=prod;                                // adding the product
    if (stoi>0) then o:='+'+o;              // insert "+" in case of positive stoi coeff
    insert(' ',o,2);                        // insert space between stoi coeff and its sign
    o:=' '+o;                               // cosmetic indent
    format_prod:=o;
    end;
end;


// - tagit body

var i, j, k, l, e, ped, d, m, d_ib : integer;
    a, b : ansistring;
    qae : array[1..2] of integer;           // q-ty of atoms in educts
    qap : integer;                          // --"-- in a given product
    upl, af, pf, pf_ib, savstoi : real;
    labbr : string;                         // current eq. abbr (of eqs[eqnfile[l].eqno].abbr)
    icode : boolean = false;                // is code added flag
{$IFDEF TAG_EXPL}
    x : integer;
{$ENDIF}

// main part

begin

write('tagit(',spcfname,':',eqnfname,',',fname,'): ');

write('/',{$IFNDEF TAG_EXPL}'NOT ',{$ENDIF}'EXPLICIT/ ');

assign(f,fname);
rewrite(f);

// some info
// writeln(f,eqnfile[1].line);
// writeln(f);
writeln(f,'#INLINE F90_GLOBAL');
writeln(f,'! -------------------------------------------------------------------------');
writeln(f,'! current mechanism (',paramstr(1),') is (isotopically) tagged by [imtag]');
writeln(f,'! configuration: ',tagname);
  write(f,'! element: <',isoelem,'>, ',_isos,' class(es)/isotopologue(s):');
for i:=1 to _isos do
    write(f,' ',clsname[i]);
writeln(f);
writeln(f,'! # of tagged / added species: ',_utsl,' / ',(_utsl+1)*_isos+1);
writeln(f,'! # of reactions in the selected mechanism: ',_eqs);
writeln(f,'! # of tagged reactions: ',nooftagreac,' (',_src,' subs)');
{$IFDEF USE_PT}
writeln(f,'! # of added passive tracers (PTs): ',_ptracs_conf,' (',_ptracs_conf-nooftagreac,' add.kie)');
{$ENDIF}
writeln(f,'! current tagging is ',{$IFNDEF TAG_EXPL}'NOT ',{$ENDIF}'EXPLICIT');
  write(f,'! composition transfer mode: ');
if (itransm=1) then writeln(f, 'single-substituted isotopologues')
               else writeln(f, 'molecular/elemental (not isotopic)');
  write(f,'! isotope exchange r-n rates are: ');
{$IFDEF IEX_ONWARD_RATES}
writeln(f,'given as onward for the first reactant');
{$ELSE}
writeln(f,'refer to the regular reaction rate');
{$ENDIF}
writeln(f,'! PTs added: ',{$IFDEF USE_PT}'YES'{$ELSE}'NO'{$ENDIF},
             ', for KIE monitoring: ',{$IFDEF USE_PT_KIE}'YES'{$ELSE}'NO'{$ENDIF},'');
writeln(f,'! =', datetimetostr(now),'=');
writeln(f,'! ----- below is the modified mechanism -----------------------------------');
writeln(f,'#ENDINLINE');
writeln(f);

for l:=1 to _eqnfile do
    if (eqnfile[l].iseq) then
       begin

       show_runner(10);  // shows that we're not stuck whilst tagging a huge eqn file

       if (eqs[eqnfile[l].eqno].itag) then           // working with equation or just passing line
          with eqs[eqnfile[l].eqno] do

          begin
          labbr:=abbr; // current eq. abbr (for inner cycles dealing with sources)
{$IFNDEF TAG_EXPL}
 {$IFNDEF IEX_REGREF}
          if (iiex) then write(f,'// ');  // disabling creation of reference reaction for isotope exchange
 {$ENDIF}
          // first, original eq with a production tracer added
          a:='<'+abbr+'> ';

          // educts
          a+=educ[1];
          if (educ[2]<>'') then a+=' + '+educ[2];
          a+=' =';

          // products
          for j:=1 to _prod do
              a+=format_prod(prod[j],stoi[j]);

 {$IFDEF USE_PT}
          if (pos(substr(ptsyntax,'@',abbr),eqnfile[l].line)=0) then
             // adding production tracer
             a+=format_prod(substr(ptsyntax,'@',abbr));
 {$ENDIF}

 {$IFDEF USE_PT_UPL}
          // checking unaccounted loss/production only in regular r-n
//        writeln(stderr,'> ',abbr,': ',qatm_educ,' -> ',qatm_prod,' = ',qatm_prod-qatm_educ);
          upl:=(qatm_prod-qatm_educ);   // qatm_* is inferred in imcom while reading eqn
          if (upl>+eps_massbal) then a+=format_prod(substr(ptpsyntax,'@','U'+cfgname+isoelem),upl);      // prod
          if (upl<-eps_massbal) then a+=format_prod(substr(ptlsyntax,'@','U'+cfgname+isoelem),abs(upl)); // loss
 {$ENDIF}

          // adding reaction rate
 {$IFNDEF USE_DKRATE}
          a+=' : '+trim(phys));
 {$ELSE}
          // retain the original markers/citation no. in the regular eqn.
          a+=' :'+' {%'+imcom_ext4marks(phys,'{%','}')+'} '
                       +substr(dkratesyntax,'@',abbr)+'; '
                       +' {&'+imcom_ext4marks(phys,'{&','}')+'}';
 {$ENDIF}

{$ENDIF}
          // adding "ever tagged" number
          a+=' {'+_etagspsy+inttostr(etag)+'}';

          // writing up
          writeln(f,substr(a,'= +','='));

          // now the tagged reactions come

          // ok, is it a isotope exchange/subsitution or normal equation for tagging?
          if (not(isrc) and not(iiex)) then
             begin

//{$ELSE}
//`             // if there is only one tagged educt and it is #2 then swapping educts in explicit tagging
//`             if (in_tsl(educ[2]) and not(in_tsl(educ[1]))) then
//`                begin a:=educ[1]; educ[1]:=educ[2]; educ[2]:=a; end;
//{$ENDIF}

             // q-ty of atoms in educts
             for e:=1 to 2 do
                 if (in_tsl(educ[e])) then qae[e]:=tsl[no_tsl(educ[e])].qatm else qae[e]:=0;

             for e:=1 to 2 do          // cycle through educts
              // current educt equation
              if in_tsl(educ[e]) then
                 for k:=1 to _isos do  // cycling through isotopologues
{$IFDEF TAG_EXPL}
                  for x:=1 to _isos do // another cycle for explicit tagging
{$ENDIF}
                    begin

                    a:=' <';
                    a+=substr(trsyntax,'@',abbr)+clsname[k];
                    if (if2t) then                       // in case both educts are tagged
{$IFNDEF TAG_EXPL}
                       if (educ[1]<>educ[2]) then
                          a+='i'+educ[e]                 // impl. quad. nomenclature: G4110I12CiSRC
                       else
                           a+='i2'+educ[e];              // impl. quad. nom. for same educts: G4110I12Ci2SRC
{$ELSE}
                       a+='e'+inttostr(x);               // expl. quad. nomenclature: G4215aI13Ce2
{$ENDIF}
                    a+='> ';

                    a+=giveiso(educ[e],k);               // isotopologue
                    if (educ[3-e]<>'') then              // '3-e' gives 2(e=1), 1(e=2)
{$IFNDEF TAG_EXPL}
                       a+=' + '+educ[3-e];               // 2nd educt,  impl.: x=0 => giveiso returns "regular" spec
{$ELSE}
                       a+=' + '+giveiso(educ[3-e],x);    // explicit: corresp. other educt isotopologue
{$ENDIF}
                    a+=' = ';                            // eq sign
{$IFNDEF TAG_EXPL}
                    a+=educ[3-e];                        // replicating "regular" @right side
{$ENDIF}

                    // accounting the atom fraction of current educt: account in stoi`s of the products
                    af:=qae[e]/(qae[1]+qae[2]);          // this is a place for improvement (now using source subs)!

                    for j:=1 to _prod do     // cycling products
                        begin
                        savstoi:=stoi[j];
{$IFNDEF TAG_EXPL}
                        if (in_tsl(prod[j])) then        // filtering out all non-tagged species in implicit
{$ENDIF}
                           begin
                         //if (a[length(a)-1]<>'=') then a+=' + '; // if right side is not empty   // switched off due to using format_prod()
{$IFNDEF TAG_EXPL}
                           stoi[j]:=stoi[j]*af;          // quadruplicated equations
{$ELSE}
                       ???       stoi[j]:=stoi[j]/_isos; // _isos-replicated equations (explicit tagging case)
{$ENDIF}
                           // calculating the transfer of the rare isotope from educt molecule to current product molecule
                           if (in_tsl(prod[j]) and (k>1) and (itransm=1)) then // if reaction is for the rare (singly subs.) isotopologue transfer
                              begin
                              qap:=tsl[no_tsl(prod[j])].qatm;          // qty. of atoms in product
                              pf:=qap/qae[e];                          // prob = A(prod)/A(educ), i.e. 1/3 for C3H6 -> CH3O2
                              // if # of atoms in product = # atoms in educt: straightforwardly creating rare isotopologue
                              // if # of atoms in product > # atoms in educt: adding atoms from abundant isotopologue pool
                              // if # of atoms in product < # atoms in educt: freeing abundant atoms excess to the abundant isotopologue pool
                              end
                           else
                               pf:=1;  // fractions, mol. tagging or only abun. isotopologue cases

                           a+=format_prod(giveiso(prod[j],k,abbr),stoi[j]*pf);       // rare isotopologue
                           a+=format_prod(giveiso(prod[j],1,abbr),stoi[j]*(1.0-pf)); // abun isotopologue, not created if pf=1
                           end;

                        // adding production PT
                        if (in_busl(prod[j])) then
                           if (in_tsl(prod[j])) then a+=format_prod(substr(ptpsyntax,'@',giveiso(prod[j],k,abbr)),stoi[j]*pf)
                                                       +format_prod(substr(ptpsyntax,'@',giveiso(prod[j],1,abbr)),stoi[j]*(1.0-pf)) // mz_sg_160209 fix: PTs when shifts are enabled
                                                else a+=format_prod(substr(ptpsyntax,'@',clsname[k]+prod[j]),stoi[j]);

                        // reverting changes to the stoichiom. coeff
                        stoi[j]:=savstoi;
                        end;

                    // loss PTs
                    // mz_sg_150809 fix: loss PTs were created in all reactions, irresp. from tagged reacting educt
                    if (if2s) then d:=2 else d:=1;
                    for i:=1 to 2 do
                        if (in_busl(educ[i])) and (not(in_tsl(educ[i])) or (i=e)) then
                           a+=format_prod(substr(ptlsyntax,'@',clsname[k]+educ[i]),d);

 {$IFDEF USE_PT_UPL}
                    // unaccounted loss/production in tagged r-ns
//                  writeln(stderr,'>> ',abbr,': ',qatm_educ,' -> ',qatm_prod,' = ',qatm_prod-qatm_educ,' af:',af);
                    upl:=af*(qatm_prod-qatm_educ);
                    if (abs(upl)>eps_massbal) then
                       begin
                       if (upl>0) then b:=ptpsyntax;
                       if (upl<0) then b:=ptlsyntax; upl:=abs(upl);
                       if ((itransm=1) and (k>1)) then   // if reaction is for the rare (singly subs.) isotopologue transfer
                          pf:=1/qae[e]                   // unaccounted loss/prod PTs always have one atom
                       else
                           pf:=1;
                       a+=format_prod(substr(b,'@','U'+cfgname+clsname[k]),upl*pf);       // rare isotopologue
                       a+=format_prod(substr(b,'@','U'+cfgname+clsname[1]),upl*(1.0-pf)); // abun isotopologue, not created if pf=1
                       end;
 {$ENDIF}

                    // checking if there are no produts on the right side (USE_PT off & destruction to nothing - quite possible)
                    b:=trim(a); if (b[length(b)]='=') then a:=a + dummyspc + ' ';
{$IFDEF USE_PT}

 {$IFDEF TAG_EXPL}
                    // PTs in explicit case: all equations have the same PT - total production
                    if (pos(substr(ptsyntax,'@',abbr),eqnfile[l].line)=0) then
                       // adding production tracer
                       a+=format_prod(substr(ptsyntax,'@',abbr),1);
 {$ENDIF}

 {$IFDEF USE_PT_KIE}
                    // add passive tracer to the the reaction which has kie
                    for j:=1 to _kie do
                        if (kie[j].imec) and (abbr=kie[j].abbr) then
                           begin
                           a+=' + '+substr(ptsyntax,'@',abbr+clsname[k]);
  {$IFNDEF TAG_EXPL}
                           if (if2t) then a+='i'+inttostr(e);
  {$ELSE}
                           if (if2t) then a+='e'+inttostr(x);
  {$ENDIF}
                           break;
                           end;
 {$ENDIF}
{$ENDIF}

                    // adding reaction rate
{$IFNDEF USE_DKRATE}
                    a+=' :'+phys;
{$ELSE}
                    a+=' :'+' {%tag:'+cfgname+'} '+substr(dkratesyntax,'@',abbr)+'; ';
{$ENDIF}

                    // in case educts are equal, skip duplicating reaction but tagging the rate
                    if (if2s) then
                       if (e=1) then
                          insert('*(2.0)',a,pos(';',a))
                       else
                           continue;

                    // checking whether the reaction has kie
                    for j:=1 to _kie do
                        if (kie[j].imec) and (abbr=kie[j].abbr) then
                           if (kie[j].isot=tsl[no_tsl(educ[e])].isos[k]) then
                              begin
                              insert('( ',a,pos('}',a)+2);            // adding left brace
                              insert(' )'+kie[j].expr,a,pos(';',a));
                              // adding right brace and expr. to the end of phys line before ;
                              end;

                    writeln(f,trimright(a));
                    end;

             end
          else   // for 'if not(isrc) and not(iiex)'
              if (isrc) then
                 // if products have special transfer/source specification
                 with src[eqs[eqnfile[l].eqno].nsrc] do
                      for e:=1 to _trans do                 // e cycles through source specification
                          if (trans[e].ib=0) then           // skipping IB transfers here
                             begin

                             // checking whether source specification from species is in the used-TSL
                             if not(in_utsl(trans[e].src)) then
                                begin
                                writeln(' <!> warning: source specification record (',tagname,': #',e,
                                        ') for <',labbr,'> refers to species not in used-TSL (',trans[e].src,')');
                                continue;
                                end;

                             // transfer info string
                             a:='// <'+labbr+'> '+trans[e].src+' = ';
                             for j:=1 to trans[e]._dst do
                                 begin
                                 if not(trans[e].dst[j].stoi=1.0) then
                                    a+=formatfloat('##0.###',trans[e].dst[j].stoi)+' ';
                                 a+=trans[e].dst[j].spec;
                                 if (j<trans[e]._dst) then
                                    a+=' + ';
                                 end;
                             a+=' : ; {%isotrans:'+isoelem+'}';
                             writeln(f,a);

                             // if a specified source is one of the educts given in the eq, creating normal r-n
                             // else creating futile destruction r-n
                             ped:=0;                // no of present educt matching src specification
                             for j:=1 to 2 do
                                 if in_tsl(educ[j]) then
                                    if (educ[j]=trans[e].src) then
                                       begin
                                       ped:=j; break;
                                       end;

                             // composing tag reactions
                             for k:=1 to _isos do
                                 begin

                                 // educts
                                 a:=' <'+substr(trsyntax,'@',labbr)+clsname[k]+'s'+trans[e].src+'> ';
                                 if (ped>0) then j:=ped else j:=1;
                                 if (educ[3-j]<>'') then
                                    a+=educ[3-j]+' + ';
                                 if (ped>0) then a+=giveiso(educ[ped],k)
                                            else a+=educ[j];
                                 a+=' = ';
                                 if (educ[3-j]<>'') then
                                    a+=educ[3-j];
                                 if (ped=0) then a+=' + '+educ[j];

                                 // resetting atomic mass balance counters
                                 if (ped>0) then qatm_educ:=tsl[no_tsl(educ[j])].qatm else qatm_educ:=0;
                                 qatm_prod:=0;

                                 // products

                                 // q-ty of atoms in educts
                                 // source specification: contibution of only one educt
                                 qae[1]:=tsl[no_tsl(trans[e].src)].qatm;
                                 qae[2]:=0;

                                 i:=length(a);    // storing the length of eq to determine if any product was written or not
                                 for j:=1 to _prod do
                                     begin
                                     savstoi:=stoi[j];
{$IFNDEF TAG_EXPL}
                                     if (in_tsl(prod[j])) then
{$ENDIF}
                                        begin

                                        // is the product in the current transfer (e) record?
                                        d:=no_trans_dst(eqs[eqnfile[l].eqno].nsrc,e,prod[j]);
                                        if (d>0) then
                                           begin

                                           // determining if there is an IB record for this particular isotopologue
                                           d_ib:=0;  // flags that IB transfer record is found
                                           for m:=1 to _trans do
                                             if (trans[m].ib=k) then      // is an IB transfer record for isotopologue (class) k
                                              begin
                                              if ( m = no_src_trans(eqs[eqnfile[l].eqno].nsrc,giveiso(trans[e].src,k)) ) then // IB transfer record matches source isotopologue
                                                 d_ib:=no_trans_dst(eqs[eqnfile[l].eqno].nsrc,m,prod[j]);                     // matching IB destination
                                              if ( d_ib>0 ) then break;   // found IB for this src->dst
                                              end;

                                           // accounting for the element transfer fraction for this product
                                           stoi[j]:=stoi[j]*trans[e].dst[d].stoi*trans[e].weight;

                                           // calculating the transfer of the rare isotope from educt molecule to current product molecule
                                           qap:=tsl[no_tsl(prod[j])].qatm;               // qty. of atoms in product
                                           if (in_tsl(prod[j]) and (k<>1) and (itransm=1)) then // if reaction is for the rare (singly subs.) isotopologue transfer
                                              if ( d_ib>0 ) then
                                                 pf:=trans[m].dst[d_ib].stoi*trans[d_ib].weight/stoi[j]    // taking the IB factor directly
                                              else
                                                  pf:=qap/qae[1]                                // prob = A(prod)/A(educ), i.e. 1/3 for C3H6 -> CH3O2
                                           else
                                               pf:=1.0;  // fractions, mol. tagging or only abun. isotopologue cases

                                           a+=format_prod(giveiso(prod[j],k,labbr),stoi[j]*pf);       // rare isotopologue
                                           a+=format_prod(giveiso(prod[j],1,labbr),stoi[j]*(1.0-pf)); // abun isotopologue, not created if pf=1

                                           qatm_prod+=stoi[j]*qap;
                                           end;
                                        end;

                                     // adding production PTs
                                     if (in_busl(prod[j])) then
                                        if not(in_tsl(prod[j])) then
                                           a+=format_prod(substr(ptpsyntax,'@',clsname[k]+prod[j]),stoi[j])
                                        else
                                            if (d>0) then
                                               a+=format_prod(substr(ptpsyntax,'@',giveiso(prod[j],k,labbr)),stoi[j]*pf)
                                                 +format_prod(substr(ptpsyntax,'@',giveiso(prod[j],1,labbr)),stoi[j]*(1.0-pf)); // mz_sg_160209 fix: PTs when shifts are enabled

                                     // reverting changes to the stoichiom. coeff
                                     stoi[j]:=savstoi;

                                     end;

                                 if (i=length(a)) then   // no products were written, i.e. src -> empty
                                    a+=' + '+dummyspc;   // accounting for destruction anyway (before was 'continue');

                                 // loss PTs
                                 if (ped>0) then
                                    begin
                                    // mz_sg_150809 fix: loss PTs were created in all reactions, irresp. from tagged reacting educt
                                    if (if2s) then d:=2 else d:=1;
                                    if (in_busl(educ[ped])) then
                                       a+=format_prod(substr(ptlsyntax,'@',clsname[k]+educ[ped]),d);
                                    end;

{$IFDEF USE_PT_UPL}
                                 // unaccounted loss/production in tagged r-ns with src
                                 // writeln(stderr,'>> ',labbr,': ',qatm_educ,' -> ',qatm_prod,' = ',qatm_prod-qatm_educ,' af:',af);
                                 upl:=(qatm_prod-qatm_educ);
                                 if (abs(upl)>eps_massbal) then
                                    begin
                                    if (upl>0) then b:=ptpsyntax;
                                    if (upl<0) then b:=ptlsyntax; upl:=abs(upl);
                                    if ((itransm=1) and (k>1)) then   // if reaction is for the rare (singly subs.) isotopologue transfer
                                       pf:=1/qae[1]                   // unaccounted loss/prod PTs always have one atom
                                    else
                                        pf:=1;
                                    a+=format_prod(substr(b,'@','U'+cfgname+clsname[k]),upl*pf);       // rare isotopologue
                                    a+=format_prod(substr(b,'@','U'+cfgname+clsname[1]),upl*(1.0-pf)); // abun isotopologue, not created if pf=1
                                    end;
{$ENDIF}


{$IFDEF USE_PT}
 {$IFDEF TAG_EXPL}
                                 // PTs in explicit case
                                 if (pos(substr(ptsyntax,'@',labbr),eqnfile[l].line)=0) then
                                    // adding production tracer
                                    a+=format_prod(substr(ptsyntax,'@',labbr),1);
 {$ENDIF}
 {$IFDEF USE_PT_KIE}
                                 // add passive tracer to the the reaction which has kie
                                 for j:=1 to _kie do
                                     if ((kie[j].imec) and (labbr=kie[j].abbr)) then
                                        begin
                                        a+=format_prod(substr(ptsyntax,'@',labbr+clsname[k],1);
                                        break;
                                        end;
 {$ENDIF}
{$ENDIF}
                                 // adding reaction rate
{$IFNDEF USE_DKRATE}
                                 a+=' :'+phys;
{$ELSE}
                                 a+=' :'+' {%tag:'+cfgname+'} '+substr(dkratesyntax,'@',labbr)+'; ';
{$ENDIF}
                                 // weighting reaction rate with corresponding isotopologue fraction
                                 //   if the source is not the one from the left side
                                 if (ped=0) then
                                    begin
                                    insert('( ',a,pos('}',a)+2); // adding left brace
                                    insert(' )*'+substr(dsrcsubsyntax,'@',clsname[k]+trans[e].src),a,pos(';',a));
                                    end;

                                 // in case educts are equal, skip duplicating reaction but doubling the rate
                                 if ((ped>0) and (educ[1]=educ[2])) then
                                    insert('*(2.0)',a,pos(';',a));

                                 // checking whether the reaction has kie
                                 for j:=1 to _kie do
                                     if (labbr=kie[j].abbr) then
                                        if (kie[j].isot=tsl[no_tsl(trans[e].src)].isos[k]) then
                                           begin
                                           insert('( ',a,pos('}',a)+2); // adding left brace
                                           insert(' )'+kie[j].expr,a,pos(';',a));
                                           // adding right brace and expr. to the end of phys line before ;
                                           end;

                                 writeln(f,trimright(a));

                                 end;

                             end
                          else   // 'if (trans[e]._dst>0)'
              else   // 'if (not(isrc)'
                  if (iiex) then
                     begin

                     // isotope exchange reaction

                     // educts
                     for e:=1 to 2 do
                       if (_isos>1) then // multiple istopologues
                         for k:=2 to _isos do
                           with iex[niex] do
                             begin

                             // new r-n abbr
                             a:=' <'+substr(trsyntax,'@',abbr)+tsl[exspec[e]].isos[k]+'> ';

                             a+=tsl[exspec[e]].isos[k]  +' + ';    // abundant
                             a+=tsl[exspec[3-e]].isos[1]+' = ';    // rare

                             // exchanging one rare isotope atom
                             a+=tsl[exspec[e]].isos[1]+' + ';      // now rare
                             a+=tsl[exspec[3-e]].isos[k];          // now abundant

                             // rate
{$IFNDEF USE_DKRATE}
                             a+=' :'+phys;
{$ELSE}
                             a+=' :'+' {%tag:'+cfgname+'} '+substr(dkratesyntax,'@',abbr)+'; ';
{$ENDIF}
                             // checking whether the reaction has kie
                             for j:=1 to _kie do
                                 if (kie[j].imec) and (abbr=kie[j].abbr) then
                                    if (kie[j].isot=tsl[exspec[e]].isos[k]) then
                                       begin
                                       insert('( ',a,pos('}',a)+2);            // adding left brace
                                       insert(' )'+kie[j].expr,a,pos(';',a));
                                       // adding right brace and expr. to the end of phys line before ;
                                       end;

{$IFNDEF IEX_ONWARD_RATES}
                             // calculating transfer probability (as nu = q(3-e) / (q(1)+q(2)))
                             // referred to the regular r-n rate, for both reactants, both directions
                             j:=tsl[exspec[1]].qatm+tsl[exspec[2]].qatm; // total on the left side
                             insert('('+inttostr(tsl[exspec[3-e]].qatm)+'.0/'+inttostr(j)+'.0)*',a,pos('}',a)+2);
{$ELSE}
                             // calculating transfer probability (as nu = q(3-e) / q(e))
                             // referred to the onward r-n rate, reverse direction
                             if (e=2) then
                                insert('('+inttostr(tsl[exspec[3-e]].qatm)+'.0/'+inttostr(tsl[exspec[e]].qatm)+'.0)*',a,pos('}',a)+2);
{$ENDIF}
                             // voila!
                             writeln(f,trimright(a));

                             end
                       else // "fractional" isotopologues (i.e. one)
                           with iex[niex] do
                             begin

                             // new r-n abbr
                             a:=' <'+substr(trsyntax,'@',abbr)+tsl[exspec[e]].isos[k]+'> ';
                             a+=tsl[exspec[e]].isos[1]  +' + ';  // rare disappeared
                             a+=tsl[exspec[3-e]].spec+' = ';     // regular (rare+rest)

                             // exchanging one rare isotope atom
                             a+=tsl[exspec[3-e]].spec+' + ';     // replicating regular (rare+rest)
                             a+=tsl[exspec[3-e]].isos[1];        // rare appeared

                             // rate
                             //   add. modified acc. to the fraction of the abundant
{$IFNDEF USE_DKRATE}
                             a+=' :'+phys;
{$ELSE}
                             a+=' :'+' {%tag:'+cfgname+'} '+substr(dkratesyntax,'@',abbr)+
{$IFDEF PRECISE_FISO_IEX}
                                     '*(1.0 - '+substr(dsrcsubsyntax,'@',clsname[1]+tsl[exspec[3-e]].spec)+')'+
{$ENDIF}
                                     '; ';
{$ENDIF}
                             // checking whether the reaction has kie
                             for j:=1 to _kie do
                                 if (kie[j].imec) and (abbr=kie[j].abbr) then
                                    if (kie[j].isot=tsl[exspec[e]].isos[k]) then
                                       begin
                                       insert('( ',a,pos('}',a)+2);            // adding left brace
                                       insert(' )'+kie[j].expr,a,pos(';',a));
                                       // adding right brace and expr. to the end of phys line before ;
                                       end;

{$IFNDEF IEX_ONWARD_RATES}
                             // calculating transfer probability (as nu = q(3-e) / (q(1)+q(2)))
                             // referred to the regular r-n rate, for both reactants, both directions
                             j:=tsl[exspec[1]].qatm+tsl[exspec[2]].qatm; // total on the left side
                             insert('('+inttostr(tsl[exspec[3-e]].qatm)+'.0/'+inttostr(j)+'.0)*',a,pos('}',a)+2);
{$ELSE}
                             // calculating transfer probability (as nu = q(3-e) / q(e))
                             // referred to the onward r-n rate, reverse direction
                             if (e=2) then
                                insert('('+inttostr(tsl[exspec[3-e]].qatm)+'.0/'+inttostr(tsl[exspec[e]].qatm)+'.0)*',a,pos('}',a)+2);
{$ENDIF}

                             // voila!
                             writeln(f,trimright(a));
                             end;

                     end;         // if (iiex)

          end
       else // 'if (...itag)'
           writeln(f,eqnfile[l].line);        // just output the equation which is not tagged

       // continuing within 'if (...iseq)' block

       // writing reactions to the end of the main block

       if (eqnfile[l].eqno=_eqs) then
          begin
{$IFDEF USE_PT_UPL}
          // LAST dummy reaction here for total loss PTs (PTL/PU*) so MECCA won't throw it away
          // et la reaction futile
          writeln(f);
          writeln(f,'{------ [imtag] - total ',isoelem,' unaccounted loss PTs dummy reaction ----------------------------------}');

          write(f,' <'+substr(trsyntax,'@',substr(ptsyntax,'@','U'+cfgname+isoelem))+'> ');
          write(f,substr(ptlsyntax,'@','U'+cfgname+isoelem),' = ',substr(ptpsyntax,'@','U'+cfgname+isoelem));
          for j:=1 to _isos do
              write(f,' + '+substr(ptpsyntax,'@','U'+cfgname+clsname[j])+' + '+substr(ptlsyntax,'@','U'+cfgname+clsname[j]));
          writeln(f,' : {%StTrG}  0.0;  {&&}');
{$ENDIF}
{$IFDEF TAG_EXPL}
          // les reactions futiles pour les substances ordinaires
          // here assumed that there are at least 3 three tagged species in the mech
          writeln(f);
          writeln(f,'{------ [imtag] - regular species dummy reactions -------------------------}');
          write('  TAG_EXPL dummy reactions: ');
          for k:=0 to ((_utsl-1) div 10) do
              begin
              write(k,' ');
              i:=  (k)*10+1;
              j:=(k+1)*10;
              if (j>_utsl) then
                 begin
                 // dec(i,j-_utsl);
                 dec(j,j-_utsl);
                 i:=j-(j mod 10)+1; if ((i+2)>_utsl) then i:=_utsl-2;
                 if (j<1) or (i<1) then
                    begin
                    writeln(' problem creating dummy reactions for regular species for explicit case. stop.');
                    halt(1);
                    end;
                 end;
              write(f,'<G000',k,'E> ',tsl[utsl[i]].spec,' + ',tsl[utsl[i+1]].spec,
                                                        ' = ',tsl[utsl[i+2]].spec);
              for x:=i+3 to j do
                  write(f,' + ',tsl[utsl[x]].spec);
              writeln(f,' : {%StTrG}  0.0;  {&&}');
              end;
          writeln;
{$ENDIF}
          end
       end    // 'if (eqnfile[l].iseq)' block
else
    begin
    // adding some code before equations
    if ((pos('#EQUATIONS',eqnfile[l].line)>0) and not(icode)) then
       begin
       icode:=true;                        // prevents adding code more than once if multiple #EQUATIONS directives are present
       if (_src>0) then                    // if source specification is present
          addinline4subs;

       if (kieproc<>_nonestr) then         // if there is KIE processing file specified
          begin
          writeln(f,'{------ [imtag] - ',tagname,' KIE processing section ----------------------}');
          write(f,imcom_parse_proc(kieproc));   // kie processing for tagging
          writeln(f,'{------ [imtag] -----------------------------------------------------------}');
          end;

// op_pj_20171021+: moved after kieproc section; otherwise KIE uninitialised variables may be used in optimised rate calculations
{$IFDEF USE_DKRATE}
       addinline4tagrates;                 // if rates calculation optimization is on (see USE_DKRATE CD)
{$ENDIF}
// op_pj_20171021-

       end;

       writeln(f,eqnfile[l].line);
    end;

writeln(f);
writeln(f);

close(f);

writeln('done');

end;

// -----------------------------------------------------------------

// additional species list - interconfigurational
// synx: :I12CO=C12 + O>
var addspecs_intr : ansistring;
   _addspecs_intr : integer;

procedure imtag_update_addspecs;
var i, j : integer;
begin
for i:=1 to _utsl do
    with tsl[utsl[i]] do
         for j:=1 to _isos do
             begin
             addspecs_intr+=':'+clsname[j]+spec+'>';
             inc(_addspecs_intr);
             end;
addspecs_intr+=':'+substr(ptlsyntax,'@',cfgname+isoelem)+'>';
for i:=1 to _isos do
    addspecs_intr+=':'+substr(ptlsyntax,'@',cfgname+clsname[i])+'>';
inc(_addspecs_intr,1+_isos);

writeln('imtag_update_addspecs: done [total ',_addspecs_intr,' specs]');
end;


// additional species list
procedure imtag_produce_speciesfile(fname : string);

var f : textfile;
    i, j, k : word;
    declspecs, a : ansistring;

// add species checking for existing in the present species list
procedure safeaddspec(name, data : string);
begin
if (in_spc(name)) then declspecs+=name+' '
                  else writeln(f,'  ',name,' = ',data);
end;


begin

write('produce_speciesfile(',fname,'): ');
declspecs:='';

assign(f,fname);
rewrite(f);

// filling in previous species file
for i:=1 to _spcfile do
    writeln(f,spcfile[i].line);

writeln(f);
writeln(f,'{-------------- [imtag] ---------------------------------------------------}');
writeln(f);
writeln(f,'{ additional species list for tagged mechanism based on: ',paramstr(1),' }');
writeln(f,'{ configuration: ',cfgname,' }');
writeln(f,'{ intended # of species to add: ',_addspecs_intr,' }');
{$IFDEF USE_PT}
writeln(f,'{ intended # of PTs to add: ',_ptracs_intr,' }');
{$ENDIF}
writeln(f);
writeln(f,'{ =',datetimetostr(now),'= }');
writeln(f);

writeln(f,'{- new atoms -----------------------------------------------------------------}');
writeln(f);
writeln(f,'#ATOMS');
writeln(f);
if (isoelem<>'') then a:=' of element '+isoelem else a:='';
for i:=1 to _isos do
    if (isomass[i]>0) then
       writeln(f,'  ',clsname[i]+isoelem,';     { mass ',isomass[i]:0:7,a,' }')
    else
       writeln(f,'  ',clsname[i]+isoelem,';     { species class ',clsname[i],'(',i,') ',a,' }');
writeln(f);

writeln(f,'{- tagged species ------------------------------------------------------------}');
writeln(f);

a:='';

(* this method is no longer used since the introduction of spc[] and work with species file
tmp:=addspecs_intr;
while (length(tmp)>0) do
      begin
      // one addspecs entry is like :C5H8>

      a:=copy2symbdel(tmp,':');   // rem ':'
      // delete(tmp,1,1); // uncomment this for fp package earlier than 5.2.0

      a:=copy2symbdel(tmp,'>');   // copy to '>'
      // delete(tmp,1,1); // uncomment this for fp package earlier than 5.2.0

      writeln(f,'  ',a,' = IGNORE;');
      end; *)

writeln(f,'#DEFVAR');
writeln(f);
for i:=1 to _utsl do
    with tsl[utsl[i]] do
         begin

         for j:=1 to _isos do
             if (trim(spc[nspc].comp[j])<>'') then
                safeaddspec(isos[j],spc[nspc].comp[j]+'; '+spc[nspc].capt[j]+'   { '+spec+' '+clsname[j]+isoelem+' isotopologue }')
             else
                 safeaddspec(isos[j],' IGNORE; '+spc[nspc].capt[j]+'   { '+spec+' '+clsname[j]+isoelem+' isotopologue, zero or no spec. elemental content given }');

         writeln(f);
         end;

writeln(f,'{ budgeting PTs }');
writeln(f);
for i:=1 to _busl do
    with busl[i] do
         begin

         // loss PTs
         if (iloss) then
            for j:=1 to _isos do
                safeaddspec( substr(ptlsyntax,'@',clsname[j]+spec),{spc[nspc].comp[j]}'IGNORE'+'; { '+spec+' '+clsname[j]+' loss PT }');
         // production PTs
         if (iprod) then
            for j:=1 to _isos do
                safeaddspec( substr(ptpsyntax,'@',clsname[j]+spec),{spc[nspc].comp[j]}'IGNORE'+'; { '+spec+' '+clsname[j]+' production PT }');

         writeln(f);
         end;


{$IFDEF USE_PT_UPL}
writeln(f,'{ unaccounted molecular/elemental production/loss PTs }');
writeln(f);
if (isoelem<>'') then a:=isoelem else a:='IGNORE';
safeaddspec(substr(ptpsyntax,'@','U'+cfgname+isoelem),a+'; { '+cfgname+' total '+isoelem+' unacc. prod. }');
safeaddspec(substr(ptlsyntax,'@','U'+cfgname+isoelem),a+'; { '+cfgname+' total '+isoelem+' unacc. loss  }');
for j:=1 to _isos do
    begin
    safeaddspec(substr(ptpsyntax,'@','U'+cfgname+clsname[j]),clsname[j]+isoelem+'; { '+cfgname+'total '+clsname[j]+isoelem+' unacc. prod. }');
    safeaddspec(substr(ptlsyntax,'@','U'+cfgname+clsname[j]),clsname[j]+isoelem+'; { '+cfgname+'total '+clsname[j]+isoelem+' unacc. loss }');
    end;
writeln(f);
{$ENDIF}

{$IFDEF USE_PT}
writeln(f,'{- passive production tracers ------------------------------------------------}');
writeln(f);

for i:=1 to _eqs do
    with eqs[i] do
         if (itag) then
            safeaddspec(substr(ptsyntax,'@',abbr),'IGNORE; { '+abbr+' reaction passive production tracer }');
writeln(f);

(* this method is no longer used since the introduction of spc[] and work with species file
tmp:=ptracs_intr;
while (length(tmp)>0) do
      begin
      // one ptracs entry is like :G4410>

      a:=copy2symbdel(tmp,':');   // rem ':'
      // delete(tmp,1,1); // uncomment this for fp package earlier than 5.2.0

      a:=copy2symbdel(tmp,'>');   // copy to '>'
      // delete(tmp,1,1); // uncomment this for fp package earlier than 5.2.0

      writeln(f,'  ',substr(ptsyntax,'@',a),' = IGNORE; { ',a,' reaction production }');
      end; *)

 {$IFDEF USE_PT_KIE}
// in case of changes refer to the imcom_update_ptracs section in imcom.inc
for k:=1 to _kie do
    if (kie[k].imec) then       // if KIE exist in this configuration
       with kie[k] do
            for j:=1 to _isos do
                if not( in_tsl(eqs[kie[k].eqno].educ[1]) and in_tsl(eqs[kie[k].eqno].educ[2]) ) then    // not a quadrupl. reaction
                   safeaddspec(abbr+clsname[j],'IGNORE; { '+abbr+' reaction '+isot+' KIE production tracer }')
                else
  {$IFNDEF TAG_EXPL}
                    begin            // in case both educts are tagged (quadrupled equation = quad. kie PTs)
                    safeaddspec(kie[k].abbr+clsname[j]+'i1',
                                'IGNORE; { '+abbr+' reaction '+isot+' KIE production tracers }');
                    safeaddspec(kie[k].abbr+clsname[j]+'i2',
                                'IGNORE;');
                    end;
  {$ELSE}
                    safeaddspec(kie[k].abbr+clsname[j]+'e'+inttostr(1),
                                'IGNORE; { '+abbr+' reaction '+isot+' KIE production tracers }');
                    for i:=2 to _isos do  // in case of explicit tagging (_iso-replicated equations with KIE)
                        safeaddspec(kie[k].abbr+clsname[j]+'e'+inttostr(i),'IGNORE;');
  {$ENDIF}

writeln(f);
 {$ENDIF}
{$ENDIF}

{$IFDEF ADD_DUMIND}
writeln(f,'{- indices of species not tagged but still present in the TSL -----------------}');
writeln(f);

for i:=1 to _tsl do
    if not(is_usedspec(tsl[i].spec)) then
       with tsl[i] do
            begin

            for j:=1 to _isos do
                safeaddspec(isos[j],'IGNORE; { dummy isotopologue }');
            writeln(f);

            end;
{$ENDIF}

writeln(f,'{-------------- [imtag] - end ---------------------------------------------}');
writeln(f);

close(f);

if (declspecs<>'') then write('(<!> warning, declined as already existing species: ',declspecs,') ');

writeln('done');

end;


// -----------------------------------------------------------------

// some useful code lines for boxmodel
procedure produce_imtag_code(tname, fname : string);

var f : textfile;



// main proc ----------------------------------------------------------

var t, finc : text;  // template, include file
    a, u, p, s : ansistring;
    ec : boolean;    // "empty" configuration flag

// main part
begin

write('produce_imtag_code(',fname,'): ');

// filling replacements
imcom_update_reps;

imcom_check_files_exist([tname, kieproc]);

// checking for an empty configuration
if (cfgname='') then ec:=true else ec:=false;

// output code
assign(f,fname);
rewrite(f);

// template file
assign(t,tname);
reset(t);

while not(eof(t)) do
      begin

      readln(t,a);
      imcom_make_reps(a);

      if pos('{>',a)>0 then
         begin
         // get control property and value letter    (e.g. ATOM:C or ISO:2 or SPEC:Ch3Br etc.)
         p:=upcase(trim(imcom_ext4marks(a,'{>',':')));
         u:=trim(imcom_ext4marks(a,':','}'));

         if not(imcom_condition(p,u)) then
            repeat
                  readln(t,a);
            until ( ((pos(('{<'),a)>0) and (pos(uppercase(p),uppercase(a))>0)) or (eof(t)) )
         else
             writeln(f,'! ----- imtag ----- condition in effect ('+p+':'+u+'): start -----');
         end
      else
          if pos('{<',a)>0 then
             begin 
             {do nothing}
             writeln(f,'! ----- imtag ----- condition in effect ('+p+':'+u+'): stop -----');
             end
          else
              if pos('{$TAG_INFO',a)>0 then
                 begin
              // writeln(f,'! ',eqnfile[1].line);
                 writeln(f,'!  source chem. mechanism files: ',paramstr(1),' (',_eqs,' reactions)');
                 writeln(f,'! # of tagged reactions/species: ',nooftagreac,' / ',_utsl,' of ',_tsl,' given');
                 end
              else
          if pos('{$CONF_PARAM',a)>0 then
             writeln(f,imcom_parameters)
          else
          if pos('{$CONF_LIST',a)>0 then
             if not(ec) then writeln(f,imcom_make_configslist(imcom_ext4marks(a,'[%','%]')))
             else
          else
          if pos('{$TRAC_DECL',a)>0 then
                 writeln(f,imcom_make_tracdecl(imcom_ext4marks(a,'[%','%]')))
          else
          if pos('{$TAG_SPECS',a)>0 then
             writeln(f,imcom_make_tagspecs_list(imcom_ext4marks(a,'[%','%]')))
          else
          if pos('{$x0',a)>0 then            // parameters are: classno, initexpression
             writeln(f,imcom_make_x0(imcom_ext4marks(a,'[%','%]'),
                                     imcom_ext4marks(a,'(%','%)')))   //_imtag_x0
          else
          if pos('{$f0',a)>0 then            // parameters are: classno, initexpression
             writeln(f,imcom_make_f0(imcom_ext4marks(a,'[%','%]'),
                                     imcom_ext4marks(a,'(%','%)')))   //_imtag_f0
          else
          if pos('{$RESET_PTs',a)>0 then
             write(f,imcom_make_resetPTs(imcom_ext4marks(a,'[%','%]'),
                                         imcom_ext4marks(a,'(%','%)')))
          else
          if (pos('{$ELSA',a)>0) then
             writeln(f,wrap(imcom_parse_eq_strarr_len(copy(a,pos('}',a)+1,length(a)-pos('}',a))),7,7))
          else
          if (pos('{$INCLUDE',a)>0) then
             begin
             s:=imcom_ext4marks(a,'<','>');
             if (fileexists(s)) then
                begin
                assign(finc,s);
                reset(finc);
                while not(eof(finc)) do
                      begin
                      readln(finc,s);
                      imcom_make_reps(s);
                      writeln(f,s);
                      end;
                close(finc);
                end
             else
                 writeln('  $INCLUDE <',s,'>: file not found. skipping.');
             end
          else
              writeln(f,a);
      end;

close(t);

close(f);

//writeln;
writeln('done');

end;

{$INCLUDE imdot.inc}

// main

var i, nconf : integer;
    l_eqnfname, l_spcfname : string;   // eqn&spc filenames from the last processed configuration

begin

if (paramcount<3) then
   begin
   writeln('>> MECCA kinetic meccanism (isotopic) tagging');
{$IFDEF TRACDEF_CHEMPROP}
   writeln('usage: ./imtag <spcfile>:<eqnfile> <chempropfile>:<processfile>   <tagging configuration(s) list>');
   writeln('   ex:           gas.spc:gas.eqn     chemprop.tbl:process_gas.tbl tagXX.cfg [tagXY.cfg ... ]');
   writeln('       (warning: *.tbl files will be altered, you may like to save the original ones)');
{$ELSE}
   writeln('usage: ./imtag <spcfile>:<eqnfile> <tracdeffile> <tagging configuration(s) list>');
   writeln('   ex:           gas.spc:gas.eqn    gas.tex       tagXX.cfg [tagXY.cfg ... ]');
{$ENDIF}
   halt;
   end;

writeln('[imtag]');
writeln('=',datetimetostr(now),'=');
writebreak;

{$IFDEF TAG_EXPL}

SORRY, EXPLICIT VERSION NEEDS RE-FORMULATION AND RE-CODING NOW

if (paramcount>4) then
   begin
   writeln('sorry, stop. current version is EXPLICIT, compiled with a TAG_EXPL switch.');
   writeln('explicit tagging can be performed only for one configuration, please check input parameters');
   halt(1);
   end;
{$ENDIF}


// inter-conf. initialization
imcom_init;
_addspecs_intr:=0;
addspecs_intr:='';

_conf:=paramcount-2; //writeln('_conf:',_conf);

for nconf:=1 to _conf do
    begin

    writeln;
    writeln('>> tagging configuration #',nconf,' of ',_conf,'...');
    writeln;

    //imcom_check_files_exist([paramstr(nconf+3),paramstr(1),paramstr(2)]);

    // read tagging config
    imcom_read_tag_config(paramstr(2+nconf));

    // read species ans equations files interpreting according to the loaded config
    if (nconf=1) then
       begin
       // spc and eqn filesnames are given now in first parameter as <spc>:<eqn>
       l_spcfname:=extractword(1,paramstr(1),[':']); // paramstr(1);
       l_eqnfname:=extractword(2,paramstr(1),[':']); // paramstr(2);
       end
    else
        begin           // next configuration is based on previously created spc/eqn files
        l_eqnfname:=eqnfname;
        l_spcfname:=spcfname;
        end;

    imcom_read_spc(l_spcfname);
    imcom_read_eqs(l_eqnfname);

    // skipping "empty" configurations
    if (cfgname='') then
       begin
       writeln;
       writeln('empty configuration detected, skipping tagging routines to code parsing');
       // copying input spc/eqn to destination spc/eqn
       fpsystem('cp '+l_spcfname+' '+spcfname);
       fpsystem('cp '+l_eqnfname+' '+eqnfname);
       continue;
       end;

    // updating inter-conf. PTs list
    imcom_update_ptracs;

    // updating inter-conf. additional specs. list
    imtag_update_addspecs;

    // tag the meccanism
    tagit(eqnfname);                      // creating equation file w/tagged mech
    imtag_produce_speciesfile(spcfname);  // creating species file for tagged mech
    for i:=1 to _form_conf do
        produce_imtag_code(form_conf[i,1],form_conf[i,2]);

    // check for a possible duplicate reactions
    imcom_check_eqn_dupes(eqnfname,'Dummy');

    conf[nconf]:=tagname;

    // tracer def. file
{$IFDEF TRACDEF_CHEMPROP}
    imcom_make_tracprop(extractword(1,paramstr(2),[':']),cmodel+'_'+tagname+'_chemprop.tbl',true);
    imcom_make_tracprop(extractword(2,paramstr(2),[':']),cmodel+'_'+tagname+'_process.tbl',false);
{$ELSE}
    imcom_make_tracprop(paramstr(2),tracdef,false);
{$ENDIF}

    // produce configuration .dot files for graphvis showing species interrelation
    if (_dots) then
       //   species names should be sparated with spaces, sources and destinations with '>'
       imdot_produce_dot_files(cmodel+'_'+tagname,
                               copy(dots,1,pos('>',dots)-1),
                               copy(dots,pos('>',dots)+1,length(dots)));

    end;

writeln;
if (_form_intr>0) then writeln('inter-configuration formers:');
for i:=1 to _form_intr do
    produce_imtag_code(form_intr[i,1],form_intr[i,2]);

writeln;
writeln('[imtag]: done');
writeln;

end.
