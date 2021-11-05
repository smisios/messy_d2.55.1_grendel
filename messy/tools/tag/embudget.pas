
// = embudget.pas ====================================================
// mz_sg_20200715 (gromov@tonnerre.mpic.de)
//
// (extended) mechanism budgeting tool, version 2.6
//
// em-budget main source module
// requirements: free pascal compiler (fpc) ver. 3 or higher
//
// [S.Gromov, MPIC, 2007-2020]
// ===================================================================

program embudget;

// - using imtag routines --------------------------------------------

{$DEFINE EMB}              // signals imcom that routines (e.g. tracdef. creation) is done for embudget
{$DEFINE IGNORE_NOSRC}     // ignoring missing sources when reading eqn (see imcom)

{$I imcom.inc}             // tagging common unit

procedure embudget_read_config(spcfile, eqnfile, ininame : string);

var i, j : integer;
    s : ansistring;

    ini : tinifile;
    info : tstringlist;

begin

writeln('reading configuration file: ',ininame);  // a ext
writeln;

// reading tagging info file
imcom_preprocess_ini(ininame);
ini:=tinifile.create('imcom.tmp');
info:=tstringlist.create;
info.clear;

cfgname:=replacestr(ininame,'.cfg','');
tagname:=cfgname;   // name of ext. budgeting cfg, p.e. emb_CIO

// ---------------------------------------------------------
// MECCA-specific params

// species' index syntax (opt., default is ind_@)
sisyntax:=extractword(1,ini.readstring('MECCA','sisyntax','ind_@'),_delims);

// passive tracers (PTs) naming syntax (opt., default is XPT@)
ptsyntax:=extractword(1,ini.readstring('MECCA','ptsyntax','XPT@'),_delims);
ptlsyntax:=extractword(1,ini.readstring('MECCA','ptlsyntax',substr(ptsyntax,'@','L')+'@'),_delims);
ptpsyntax:=extractword(1,ini.readstring('MECCA','ptpsyntax',substr(ptsyntax,'@','P')+'@'),_delims);

// tagged reac. naming syntax (opt., default is TAG@)
trsyntax:=extractword(1,ini.readstring('MECCA','trsyntax','TAG@'),_delims);

// ---------------------------------------------------------
// reading species to budget
writeln('- species: budgeting categories');

info.clear;
ini.readsection('SPC',info);
if (info.count<=0) then
   imcom_error('error: [SPC] section information is missing. stop.');

// initialising data
fillchar(tsl,sizeof(tsl),0); fillchar(tsla,sizeof(tsla),0);

_tsl:=0;
for j:=0 to info.count-1 do
    if (info[j]<>'') then   // necessary condition to avoid empty keys of occasional trash in cfg
       begin

       // cheating imcom
       inc(_tsl);
       tsl[_tsl].spec:=extractword(1,info[j],_delims);
       tsl[_tsl].qatm:=strtointdef(extractword(2,info[j],_delims),1);

       s:=ini.readstring('SPC',info[j],'');
       if (pos('L',upcase(s))>0) then tsl[_tsl].iloss:=true;
       if (pos('P',upcase(s))>0) then tsl[_tsl].iprod:=true;

//     writeln('read: <',info[j],'> = <',s,'> / interpreted: <',tsl[_tsl].spec,'> <',tsl[_tsl].qatm,'> <',tsl[_tsl].iloss,'> <',tsl[_tsl].iprod,'>');
       end;

// ---------------------------------------------------------
// reading budgeting groups into ncat[_cats], conditions into ccat[_cats]

for i:=1 to _tsl do
 with tsl[i],tsla[i] do
    begin

    write('  ',spec,':');

    info.clear;
    ini.readsection(spec,info);
    if (info.count<=0) then
        imcom_error('species ['+spec+'] section information is missing. stop.');

    tsla[i]._cats:=0;
    for j:=0 to info.count-1 do
        if (info[j]<>'') then   // necessary condition to avoid empty keys of occasional trash in cfg
           begin

           // cheating imcom further
           inc(_cats);
           ncat[_cats]:=info[j];
           ccat[_cats]:=ini.readstring(spec,info[j],'');

           write(' ',ncat[_cats],'(',ccat[_cats],')');
           end;

    writeln;
    end;

writebreak;

// ---------------------------------------------------------
// additional species
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
   writeln('- species to add:');
   for i:=1 to _addspc do
       writeln('  ',addspc[i]);
   writeln;
   end;

// additional equations
_addeqs:=0; fillchar(addeqs,sizeof(addeqs),0);
info.clear; ini.readsection('ADD:EQN',info);
for i:=0 to info.count-1 do
    begin
    inc(_addeqs);
    s:=ini.readstring('ADD:EQN',info[i],'');
    addeqs[_addeqs]:=trim(info[i]);
    if (s<>'') then
      addeqs[_addeqs]+=' = '+s;
    end;
if (_addeqs>0) then
   begin
   writeln('- equations to add:');
   for i:=1 to _addeqs do
       writeln('  ',addeqs[i]);
   writeln;
   end;

// reading species data using imcom
imcom_read_spc(spcfile);
// reading equations using imcom
imcom_read_eqs(eqnfile);

// cleanup
info.destroy;
ini.destroy;

end;


// budgeting routine, produces eqn/spc file

procedure budgetit(spcname, eqnname : string);

// - budgetit body

var f : text;
    i, j, k, l, e, npt : integer;
    a, b, c, out : ansistring;

   _bs : integer;
    bs : array[1..max_prod] of record
         spec : nstr;    // species to budget
         stoi : real;    // coefficient
         il : boolean;   // is loss?
         end;

begin

write('budgetit(',eqnname,':',spcname,',',cfgname,'): ');

out:='';  // output

for l:=1 to _eqnfile do
    if not((eqnfile[l].iseq) and (eqs[eqnfile[l].eqno].itag)) then
       out+=eqnfile[l].line+_LF         // just output the line if it is not an equation
    else
        with eqs[eqnfile[l].eqno] do
          begin
          c:=eqnfile[l].line;

          _bs:=0;

          // processing educts
          inc(_bs); bs[_bs].spec:=educ[1]; bs[_bs].stoi:=1.0; bs[_bs].il:=true;

          if (educ[2]<>'') then
             begin

             if (educ[1]=educ[2]) then
                bs[_bs].stoi:=2
             else
                 begin
                 inc(_bs); bs[_bs].spec:=educ[2]; bs[_bs].stoi:=1.0; bs[_bs].il:=true;
                 end;
             end;

          // products
          for j:=1 to _prod do
              begin
              inc(_bs); bs[_bs].spec:=prod[j]; bs[_bs].stoi:=stoi[j]; bs[_bs].il:=false;
              end;

          // checking if we need to add PTs for species in bs[]
          a:='';
          for k:=1 to _bs do
              if in_tsl(bs[k].spec) then
                 with tsl[no_tsl(bs[k].spec)],tsla[no_tsl(bs[k].spec)] do
                      begin

                      // loss PTs
                      if (iloss) and (bs[k].il) then
                         begin

                         if (abs(bs[k].stoi)<>1.0) then b:=floattostr(abs(bs[k].stoi))+' ' else b:='';
                         if not(bs[k].stoi<0.0) then b:=' + '+b else b:=' - '+b;

                         // checking if r-n in a certain PT category
                         for e:=1 to _cats do
                          for i:=1 to wordcount(ccat[e],_delims) do
                             if (pattmatch(abbr, extractword(i,ccat[e],_delims), true)) then
                                begin
                                a+=b+substr(ptlsyntax,'@',ncat[e]+spec);   // ' + '
                                uloss[e]:=true;    // flag if used
                                end;
                         end;

                      // production PTs
                      if (iprod) and not(bs[k].il) then
                         begin

                         if (abs(bs[k].stoi)<>1.0) then b:=floattostr(abs(bs[k].stoi))+' ' else b:='';
                         if not(bs[k].stoi<0.0) then b:=' + '+b else b:=' - '+b;

                         // checking if r-n in a certain PT category
                         for e:=1 to _cats do
                          for i:=1 to wordcount(ccat[e],_delims) do
                             if (pattmatch(abbr, extractword(i,ccat[e],_delims), true)) then
                                begin
                                a+=b+substr(ptpsyntax,'@',ncat[e]+spec);   // ' + '
                                uprod[e]:=true;    // flag if sed
                                end;
                         end;
                      end;

          // inserting PTs before the reaction rate
          insert(' '+trim(a)+' ',c,pos(':',c));   // inserting in original EQ

          out+=c+_LF;
          end;

// some info
a:='';
a+='#INLINE F90_GLOBAL'+_LF;
a+='! -------------------------------------------------------------------------'+_LF;
a+='! current mechanism is budgeted by [embudget]'+_LF;
a+='! configuration: '+cfgname+_LF;
a+='! list of budgeted species/categories:'+_LF;
npt:=0;
for i:=1 to _tsl do
 with tsl[i],tsla[i] do
    begin
    a+='!  prod '+spec+' :';
    for j:=1 to _cats do
        if (uprod[j]) then
           begin
           a+=' '+ncat[j]; //+'('+ccat[j]+') ';
           inc(npt);
           end;
    a+=_LF;
    a+='!  loss '+spec+' :';
    for j:=1 to _cats do
        if (uloss[j]) then
           begin
           a+=' '+ncat[j]; //+'('+ccat[j]+') ';
           inc(npt);
           end;
    a+=_LF;
    end;
a+='! # of added passive tracers (PTs): '+inttostr(npt)+_LF;
a+='! # of reactions in the selected mechanism: '+inttostr(_eqs)+_LF;
a+='! # of budgeted reactions: '+inttostr(nooftagreac)+_LF;
a+='! ='+datetimetostr(now)+'='+_LF;
a+='! ----- below is the modified mechanism -----------------------------------'+_LF;
a+='#ENDINLINE'+_LF+_LF;

out:=a+out;

// flushing eqn
assign(f,eqnname);
rewrite(f);
write(f,out);
close(f);

// adding PTs to the SPC
assign(f,spcname);
rewrite(f);

// flushing original SPC
for i:=1 to _spcfile do
    writeln(f,spcfile[i].line);

writeln(f);
writeln(f,'{-------------- [embudget] ----------------------------------------------------}');
writeln(f);
writeln(f,'{ additional passive tracers (PTs) for mechanism: ',eqnname,' }');
writeln(f,'{ configuration: ',cfgname,' }');
writeln(f,'{ # of PTs added: ',npt,' }');
writeln(f);
writeln(f,'{ =',datetimetostr(now),'= }');
writeln(f);
writeln(f,'{- passive production/loss tracers --------------------------------------------}');
writeln(f);

for i:=1 to _tsl do
    with tsl[i],tsla[i] do
         begin
         writeln(f,'{--- ',spec,' ---}');
         for j:=1 to _cats do
             begin
             if (iloss) and (uloss[j]) then
                writeln(f,'  ',substr(ptlsyntax,'@',ncat[j]+spec),' = IGNORE; ',
                          '{ '+ncat[j]+' ('+ccat[j]+') category loss passive tracer }');
             if (iprod) and (uprod[j]) then
                writeln(f,'  ',substr(ptpsyntax,'@',ncat[j]+spec),' = IGNORE; ',
                          '{ '+ncat[j]+' ('+ccat[j]+') category production passive tracer }');
             end;
         writeln(f);
         end;

writeln(f,'{-------------- [embudget] - end ----------------------------------------------}');
writeln(f);

close(f);

writeln('done');

end;


var nconf : integer;
    l_eqnfname, l_spcfname : string;   // eqn&spc filenames from the last processed configuration

const genname = 'messy_mecca_tag_embudget';
      spcname = genname+'.spc';
      eqnname = genname+'.eqn';
      texname = genname+'.tex';

begin

imcom_init();

if (paramcount<3) then
   begin
   writeln('>> MECCA kinetic meccanism extended budgeting');
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

writeln('[embudget]');
writeln('=',datetimetostr(now),'=');
writebreak;

_conf:=paramcount-2;

for nconf:=1 to _conf do
    begin

    writeln;
    writeln('>> budgeting ',nconf,' of ',_conf,' configuration(s)...');
    writeln;

    // read species ans equations files interpreting according to the loaded config
    if (nconf=1) then
       begin
       // spc and eqn filesnames are given now in first parameter as <spc>:<eqn>
       l_spcfname:=extractword(1,paramstr(1),[':']); // paramstr(1);
       l_eqnfname:=extractword(2,paramstr(1),[':']); // paramstr(2);
       end
    else
        begin           // next configuration is based on previously created spc/eqn files
        l_eqnfname:=eqnname;
        l_spcfname:=spcname;
        end;

    embudget_read_config(l_spcfname, l_eqnfname,
                         paramstr(2+nconf));
    budgetit(spcname, eqnname);

{$IFDEF TRACDEF_CHEMPROP}
    imcom_make_tracprop(extractword(1,paramstr(2),[':']),genname+'_chemprop.tbl',true,(nconf>1));
    imcom_make_tracprop(extractword(2,paramstr(2),[':']),genname+'_process.tbl',false,(nconf>1));
{$ELSE}
    imcom_make_tracprop(paramstr(2),texname,false,(nconf>1));
{$ENDIF}

    end;

writeln;
writeln('[embudget]: done');
writeln;

end.

