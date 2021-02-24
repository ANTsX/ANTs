---
name: Bugs or other run time problems
about: Report potential bugs or other problems running ANTs programs.

---

**Before opening an issue**

Please check the ANTs wiki:

https://github.com/ANTsX/ANTs/wiki

The Github wiki search only covers page titles, but you can do a full 
text search by entering into Google:

```
<your search terms> site:https://github.com/ANTsX/ANTs/wiki
```

Please check previous issues on Github at

https://github.com/ANTsX/ANTs/issues

```
is:issue <your search terms>
```

If you need to open an issue, you can remove this text and fill out the sections
below. Please leave the section titles in place, but replace the instruction
text with details of your issue.


**Describe the problem**
Please describe the problem.


**To Reproduce**
Steps to reproduce the problem. The majority of issues with ANTs are specific to
the data, so uploading example data will make it much easier to provide help. If
it is impossible to share the data in question, attempting to reproduce the
problem with other public data is helpful.

Please include:

 * The exact command line as it was run. Please run all commands with
   verbose output where possible. If your command takes a long time, please try
   to produce a faster example that still shows the problem (eg, by running
   fewer iterations). 

 * The full verbose output printed to the terminal when you run the command. If
   this is very long, please paste into a text file and upload as an attachment.
 
 * Example input data that goes with the command line you provide. You
   can share your own data or link to any public data that reproduces the
   problem you are experiencing. If you cannot reproduce the problem with public
   data, please let us know which data sets you tried.

If uploading data as an attachment, please try to minimize the file size by
compressing, downsampling or otherwise creating smaller images that demonstrate
the problem.


**System information (please complete the following information)**
Many issues are specific to a particular system. Please include all information
about your computing environment.

 - OS: [e.g. Mac OS]
 - OS version: [e.g. 10.15.1]
 - Type of system: [Desktop, laptop, HPC cluster, cloud instance,
   other]
 
   If you are building inside a virtual machine, container, Cygwin, 
   Windows Subsystem for Linux, or other non-native environment, please 
   let us know and include details of both the virtual Linux and the 
   host OS.

**ANTs version information**
 - ANTs code version: [output of antsRegistration --version]
 - ANTs installation type: [Compiled from source, downloaded binary
   (from where?), installed by another package (which?), installed in
   container (URL?), other (please specify)]

**Additional information**
Add any other information that might help solve the problem.
