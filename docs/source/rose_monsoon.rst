
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
