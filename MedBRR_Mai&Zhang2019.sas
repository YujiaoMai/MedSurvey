
/***************************************************************/;
/* SAS macro of the BRR test for mediation effects             */;
/* Authors: Yujiao Mai and Hui Zhang                           */;
/* Affiliation: St. Jude Children's Research Hospital          */;
/* Email: yujiao.mai@stjude.org and hui.zhang@stjude.org       */;
/* Licence: MIT                                                */;
/* Date: 2018-2019                                             */;
/* URL: https://github.com/YujiaoMai/MedSurvey/                 */;
/***************************************************************/;

OPTIONS PS=58 LS=80 NODATE NONUMBER;

%macro rma(ina=,RR=160);
	%do i=0 %to &RR;
		%if %sysfunc(exist( &ina&i,data)) %then %do;
		  proc sql;
		    drop table &ina&i;
		  quit;
		%end;
	%end;
%mend;
%macro rvo(di=);/*Remove DATASET*/
	%if %sysfunc(exist(&di,data)) %then %do;
	  proc sql;
	    drop table &di;
	  quit;
	%end;
%mend rvo;

**** generate the model syntax; 
%macro fitmodel(xvar=, yvar=, mvars=,zvars=,delim='/',datain=,rwname=,r=1);%macro _; %mend _;
	%rvo(di=Modelsyntax);
	Proc iml;
	mnames=scan(&mvars, 1:( 1 + countc(&mvars, &delim)));
	znames=scan(&zvars, 1:( 1 + countc(&zvars, &delim)));
	nm = prod(dimension(mnames));
	nz = prod(dimension(znames));
	tn=nm+4;
	Blanks = BlankStr(200);
	modelmatrix=J(tn,4,Blanks);	
	Procstr = BlankStr(200); Procstr=cat('PROC CALIS data=',&datain,' OUTEST=ESTS',&r,' METHOD=ML NOPRINT;');
	weightstr = BlankStr(100); weightstr=cat('WEIGHT ',&rwname,&r,'; Run;');
	vars = BlankStr(200); vars=cat('VAR ',&yvar,' ',&xvar); 
	do i=1 to nm;vars=cat( vars,' ',mnames[i]); end; do j=1 to nz;vars=cat( vars,' ',znames[j]); end; 
	vars=cat(vars,';'); 
	tn=1;
	modelmatrix[tn,1]=vars;modelmatrix[tn,2]=char(tn); 
		modelmatrix[tn,3]=Procstr;
	tn=tn+1;
	modelmatrix[tn,1]='LINEQS ';modelmatrix[tn,2]=char(tn);
	eqsy=cat(&yvar,'=','u',0,'*Intercept','+','c',0,'*',&xvar);
	do i=1 to nm;eqsy=cat(eqsy,'+','b',i,'*',mnames[i]);end;
	do j=1 to nz;eqsy=cat(eqsy,'+','d',0,j,'*',znames[j]);end;
	eqsy=cat(eqsy,'+','e',0);
	tn=tn+1;
	modelmatrix[tn,1] = cat(eqsy,',');modelmatrix[tn,2]=char(tn);
	do i=1 to nm;
		eqs = cat(mnames[i],'=','u',i,'*Intercept','+','a',i,'*',&xvar);
		do j=1 to nz; eqs = catt(eqs,'+','d',i,j,'*',znames[j]);end;
		eqs = cat(eqs,'+','e',i);
		if i=nm then eqs=cat(eqs,';');else eqs=cat(eqs,',');
		tn=tn+1;
		modelmatrix[tn,1]=eqs;modelmatrix[tn,2]=char(tn);
	end;	
	if nm>=2 then do; allc=ALLCOMB(nm, 2);ncomb=nrow(allc); end;
	else ncomb=0;	
	covs='COV ';
		do i=1 to ncomb;
			if (i=1) then do; covs=cat(covs,'e',allc[i,1],' ', 'e',allc[i,2],'=','cov',allc[i,1],allc[i,2]); end;
			else do; covs=cat(covs,', ','e',allc[i,1],' ', 'e',allc[i,2],'=','cov',allc[i,1],allc[i,2]); end;
		end;
		if covs='COV ' then do; covs='';end;
		covs=cat(covs,';');
	tn=tn+1;
	modelmatrix[tn,1]=covs;modelmatrix[tn,2]=char(tn);modelmatrix[tn,4]=weightstr;
	mattrib modelmatrix colname={'syntax' 'linenum' 'procs' 'weights'};	
	create modelsyntax from modelmatrix[colname={'syntax' 'linenum' 'procs' 'weights'}];append from modelmatrix;
	/*show contents;*//*close modelsyntax;*/
	quit;
	proc sql noprint;
		select distinct 
		   cat(procs,syntax,weights) 
		   length=5000 into:model separated by ' ' 
	  	from Modelsyntax ORDER BY linenum;
	quit;

	%macro runSteps; &model.;%mend; 
	%runSteps;
%mend;
/*%fitmodel(xvar='workban', yvar='numcg', mvars='sp_adltban/sp_kidsban',zvars='PRTAGE/NumKid',delim='/',*/
/*		datain='MedData', rwname='repwgt', r=10);*/
%macro fitmodels(xvar=, yvar=, mvars=,zvars=,delim='/',datain=,rwname=,RR=, tableout=);
		%rvo(di=&tableout);
		%DO r=0 %TO &RR;			
			%fitmodel(xvar=&xvar, yvar=&yvar, mvars=&mvars,zvars=&zvars,delim=&delim,datain=&datain,rwname=&rwname, r=&r);
			%IF &r=0 %Then %DO; DATA &tableout; SET ESTS&r; WHERE _TYPE_="PARMS"; RUN; 
			%END; %ELSE %DO; proc append base=&tableout data=ESTS&r; WHERE _TYPE_="PARMS"; RUN;%END;
		%END;
%mend;
%macro parmtests(estparm=,parnames='a0/d01/d02/d11/d12/d21/d22/v0/v1/v2',Fay=4,RR=160);%macro _; %mend _;
	Proc iml;
		use &estparm; read all var _ALL_ into estall[colname=varNames];names=varNames;close &estparm;
		nameparms=scan(&parnames, 1:( 1 + countc(&parnames, '/')));
		estparms=estall[,nameparms];
		R=&RR;Fay=&Fay;
		estbar=estparms[1,];estrps=estparms[2:R+1,]; 
		estds=estrps-estbar ; /* Y_r-Y_0 */
		temp1=(estds#estds); /* (Y_r-Y_0)^2 */ 
		sds=SQRT(FAY*MEAN(temp1)); /* Sqrt(Fay*sum((Y_r-Y_0)^2)/R) */
		zs=estbar/sds;/* Z statistics */
		ps=(1-probnorm(abs(zs)))*2;/* two-sided p-values */
		tests=estbar//sds//zs//ps;
		mattrib tests colname=nameparms rowname={'est' 'SD_BRR' 'Z' 'p-value'};
		print tests[format=8.5];
	QUIT;
%mend;
%macro medtests(estparm=,anames='a1/a2',bnames='b1/b2',mednames='sp_adltban/sp_kidban',Fay=4,RR=160,adjmethod=holm);%macro _; %mend _;
	Proc iml;
	use &estparm; read all var _ALL_ into estall[colname=varNames];names=varNames;close &estparm;
	parmas=scan(&anames, 1:( 1 + countc(&anames, '/')));	parmbs=scan(&bnames, 1:( 1 + countc(&bnames, '/'))); 
	mednames=scan(&mednames, 1:( 1 + countc(&mednames, '/')));
	as=estall[,parmas];bs=estall[,parmbs];
	nas=ncol(parmas);nbs=ncol(parmbs);
	R=&RR;Fay=&Fay;
	IF nas=nbs THEN
	/*BEGIN: Calculate the BRR standard errors and Z tests */
		abs=as#bs; /* Construct the product ab*/
		ests=as||bs||abs; 			
		estbar=ests[1,];estrps=ests[2:R+1,]; /* Separate the estimate0 (Y_0) and the esimate1-160 (Y_r)*/
		estds=estrps-estbar; /* Y_r-Y_0 */
		temp1=(estds#estds); /* (Y_r-Y_0)^2 */ /* Equivalent to temp2=(estds##2) */ 
		sds=SQRT(FAY*MEAN(temp1)); /* Sqrt(Fay*sum((Y_r-Y_0)^2)/R) */
		zs=estbar/sds; /* Z statistics */
		ps=(1-probnorm(abs(zs)))*2; /* two-sided p-values */
	/*BEGIN: Format printing*/
		idxa=1:nas; idxb=nas+1:nas+nbs; idxab=nas+nbs+1:nas+nbs+nbs;
		estas=estbar[,idxa];estbs=estbar[,idxb]; estabs=estbar[,idxab];
		sdas=sds[,idxa];sdbs=sds[,idxb];sdabs=sds[,idxab];
		zas=zs[,idxa];zbs=zs[,idxb];zabs=zs[,idxab];
		pas=ps[,idxa];pbs=ps[,idxb];pabs=ps[,idxab];		
		sdSB=SQRT((estas##2)#(sdbs##2)+(estbs##2)#(sdas##2));zSB=estabs/sdSB;pSB=(1-probnorm(abs(zSB)))*2;
		tests=estas`||sdas`||zas`||pas`||estbs`||sdbs`||zbs`||pbs`||estabs`||sdabs`||zabs`||pabs`;
		mattrib tests colname={'a' 'SD_BRR(a)' 'Z(a)' 'p-value(a)' 'b' 'SD_BRR(b)' 'Z(b)' 'p-value(b)' 
			'ab' 'SD_BRR(ab)' 'Z(ab)' 'p-value(ab)'} rowname=mednames;
/*		print tests[format=8.5];		*/
		medvars = mednames`;mattrib medvars colname={'mediator' 'test'};
		create medtable from tests;append from tests;close medtable;
		create mediator from medvars[colname={'mediator'}];append from medvars;close mediator;	
	QUIT;
	data abpvalues(rename=(COL12=RAW_P)); set medtable; keep COL12;run;
	Data medtable(rename=(COL1=a COL2=SD_BRR_a COL3=Z_a COL4=pValue_a COL5=b COL6=SD_BRR_b 
							COL7=Z_b COL8=pValue_b COL9=ab COL10=SD_BRR_ab COL11=Z_ab COL12=pValue_ab));
		set mediator;set medtable;
		label COL1="a" COL2="SD_BRR(a)" COL3="Z(a)" COL4="p-value(a)" COL5="b" COL6="SD_BRR(b)" 
							COL7="Z(b)" COL8="p-value(b)" COL9="ab" COL10="SD_BRR(ab)" COL11="Z(ab)" COL12="p-value(ab)";
	run;
	ODS LISTING CLOSE;ODS HTML CLOSE;	
	proc multtest inpvalues=abpvalues &adjmethod;ods output pValues=adjps;run; 
	data adjpValues(rename=(RAW_P=pValueRaw)); set mediator; set abpvalues; set adjps; drop Raw test;run;	
	ODS html file='medtest.html';
	title "Tests for Mediation Effects";
	proc print data=medtable;run;
	title "Adjusted p-values";
	proc print data=adjpValues;run;
	%rvo(di=Mediator);%rvo(di=abpvalues);%rvo(di=Adjps);
%mend;
%macro MediationBRR(xvar=, yvar=, mvars=,zvars=,delim='/',datain=,rwname=,RR=10,Fay=4,adjmethod=holm);%macro _; %mend _;
	ODS LISTING CLOSE;ods html close;
	%fitmodels(xvar=&xvar, yvar=&yvar, mvars=&mvars,zvars=&zvars,delim=&delim,datain=&datain,rwname=&rwname,RR=&RR, tableout=ESTSparam);
	%rma(ina=Ests,RR=&RR);
	proc iml;		
		mnames=scan(&mvars, 1:( 1 + countc(&mvars, &delim)));
		nm = prod(dimension(mnames));
		anames=BlankStr(200);bnames=BlankStr(200);medeffs = BlankStr(200);
		anames=cat("'",'a',1); bnames=cat("'",'b',1); medeffs=cat("'",'a',1,'b',1);
		if nm >=2 then do;
			do i=2 to nm; anames=cat(anames,'/','a',i); bnames=cat(bnames,'/','b',i); medeffs=cat(medeffs,'/','a',i,'b',i); end;
		end;
		anames=cat(anames,"'"); bnames=cat(bnames,"'"); medeffs=cat(medeffs,"'");
		call symputx("as", anames);            /* create macro variable */	
		call symputx("bs", bnames);
		call symputx("abs", medeffs);	
	quit;	
	%medtests(estparm=Estsparam,anames=&as,bnames=&bs,mednames=&mvars,Fay=&Fay,RR=&RR,adjmethod=&adjmethod);	
%mend;

/** Example: Data 'PisaMed' can be download from https://github.com/YujiaoMai/MedSurvey/ **/;
%MediationBRR(xvar='ParenSpt', yvar='Math', mvars='StuMtv',zvars='AGE/SEX/StuAnxt',delim='/',
		datain='PisaMed',rwname='W_FSTURWT',RR=80,Fay=4);
/** Example: Data 'TUSMed' can be download from https://github.com/YujiaoMai/MedSurvey/ **/;
%MediationBRR(xvar='workban', yvar='numcg', mvars='sp_adltban/sp_kidsban/sp_homeban',zvars='PRTAGE/PESEX/NumKid',delim='/',
		datain='TUSMed',rwname='repwgt',RR=160,Fay=4,adjmethod=HOLM HOMMEL FDR);



		

