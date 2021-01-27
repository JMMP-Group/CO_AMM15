
====================================================================
Notes on rose suites and  associated branches for running on Monsoon
====================================================================

Here we keep for the record steps to create a nominal package branch on 
the Paris svn repository for the running of AMM15 suites.

We start at 4.0.2 for now but will update once the GO package branch is chilled


Background to rose suite and package branch
==================================================

See: 
:ref:`CO9p0_AMM15` 
for background to experimental set up
that derives the COp0_AMM15 suite.

Here we take that suite and create a clean suite for use on monsoon.

That has 2 basic requirements

  * Create a version of the suite on the collaboration repository
  * Create a branch on the Paris repository that the rose suite will use as the src code 

Note later experiments that incorporate Nico's changes will use this branch and suite (with namelist settings set within)
as the reference point




Create reference version of the suite under collaboration repository
=====================================================================


   
Create package branch on Paris Repository
=====================================================================
1. Create 4.0.2 package branch, starting from G08 version of the same:
 Created Ticket for 4.0.2 starter package for c09

svn copy svn+ssh://deazer@forge.ipsl.jussieu.fr/ipsl/forge/projets/nemo/svn/NEMO/branches/UKMO/NEMO_4.0.2_GO8_package svn+ssh://deazer@forge.ipsl.jussieu.fr/ipsl/forge/projets/nemo/svn/NEMO/branches/UKMO/NEMO_4.0.2_CO9_package -m "CO9_AMM15 package branch relative to 4.0.2 GO8_package_branch"

See RMED Ticket:

https://code.metoffice.gov.uk/trac/rmed/ticket/133#ticket



Move to Version 4.0.4
=======================

Also see the issues tracker:

https://github.com/JMMP-Group/CO9_AMM15/issues/15

Have a version of the rose suite at version 4.0.2 on monsoon.

The suite is: u-ca880
It can be found on the Met office collaboration repository under:
https://code.metoffice.gov.uk/svn/roses-u/c/a/8/8/0/trunk

To check it out and run: (once you have set up collab repos credentials)

    fcm co https://code.metoffice.gov.uk/svn/roses-u/c/a/8/8/0/trunk u-ca880
    cd u-ca880
    rose suite-run

Making Version 404 should be straightforward.
Note for some reason on monsoon the initial stage of nemo takes a very long time,
this is even the case if using the same executable and xios executable from the internal MO cray.
This adds considerably to the overall run time. As the suite is set up to run in a daily cycle this is a bit of an issue for longer runs.

Either we solve what is the initialization problem on monsoon or we rewrite the suite to cycle on longer timescales e.g. months.
Otherwise it probably isn't very efficient to run long runs on monsoon.
Also have 4.0.4 version of this suite:
u-ca907

uses the 404 package branch:
nemo_sources=branches/UKMO/NEMO_4.0.4_CO9_package

As in 402 we have an odd slow down (factor of ten) for the initialisation step of NEMO compared to running on
internal cray.
Source of issue as yet not known

For longer cycle runs e.g. a month the problem wont be to bad but it is a problem for daily cycles which the suite is defaulted to
adding order 2 minutes to a <8 minute run.


Rose Suites on Monsoon
=======================

u-ca996: Tests 90s times step, monthly cycle with default 4.0.4 base starting 1990 with minimal output

u-cb676: As u-ca996, use restart from 4.0.2 (RESTART_FROM_U-BP634_BASED_TEST_NEMO4.0.2_CO7_MATCH_SAFEBUILD_TIDE_V1_GLS_EN_BUG_FIX_HALF_MASK) run from 2004 tests new outputs
