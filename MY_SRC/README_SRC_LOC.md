# Trace of SRC Code

Src code comes from the branch:

https://code.metoffice.gov.uk/svn/nemo/NEMO/branches/NERC/NEMO_4.0.4_CO9_package_tides/src

Revision: 16135


This itself was branched from a GO8_package_branch  which comes from NEMO 404
we go back to when that branch dealt with the svn keywords and take all the change since then


```
------------------------------------------------------------------------
r14075 | cguiavarch | 2020-12-04 10:02:42 +0000 (Fri, 04 Dec 2020) | 2 lines

UKMO/NEMO_4.0.4_mirror : Remove SVN keywords. 
```

That branch came from NEMO 404:

```
------------------------------------------------------------------------
r13653 | smasson | 2020-10-21 13:54:07 +0100 (Wed, 21 Oct 2020) | 1 line

create release r4.0.4
```


differences since the keywords were cleared out are:


svn diff --summarize -r14075

```
M       ICE/ice.F90
M       ICE/icerst.F90
M       ICE/icethd_pnd.F90
M       ICE/icewri.F90
 M      ICE
M       OCE/BDY/bdy_oce.F90
M       OCE/BDY/bdydta.F90
M       OCE/BDY/bdyini.F90
M       OCE/BDY/bdytides.F90
 M      OCE/BDY
M       OCE/DIA/diawri.F90
A       OCE/DIA/diaprod.F90
 M      OCE/DIA
M       OCE/DOM/dom_oce.F90
M       OCE/DOM/domain.F90
M       OCE/DOM/dommsk.F90
M       OCE/DOM/dtatsd.F90
 M      OCE/DOM
M       OCE/DYN/dynspg_ts.F90
 M      OCE/DYN
M       OCE/ICB/icbrst.F90
 M      OCE/ICB
M       OCE/IOM/in_out_manager.F90
M       OCE/IOM/iom.F90
M       OCE/IOM/restart.F90
 M      OCE/IOM
M       OCE/SBC/sbcssm.F90
M       OCE/SBC/sbctide.F90
M       OCE/SBC/tide.h90
M       OCE/SBC/tide_mod.F90
M       OCE/SBC/tideini.F90
 M      OCE/SBC
M       OCE/ZDF/zdfgls.F90
M       OCE/ZDF/zdfmxl.F90
 M      OCE/ZDF
M       OCE/step.F90
M       OCE/step_oce.F90
A       OCE/test_svn_relocate.txt
 M      OCE
M       SAS/sbcssm.F90
 M      SAS
 M      .
```

